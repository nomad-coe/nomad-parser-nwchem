from __future__ import absolute_import
from nomadcore.simple_parser import SimpleMatcher as SM
from nomadcore.caching_backend import CachingLevel
from nomadcore.baseclasses import MainHierarchicalParser, CacheService
import re
import logging
import numpy as np
LOGGER = logging.getLogger("nomad")


#===============================================================================
class NWChemMainParser(MainHierarchicalParser):
    """The main parser class that is called for all run types. Parses the NWChem
    output file.
    """
    def __init__(self, file_path, parser_context):
        """
        """
        super(NWChemMainParser, self).__init__(file_path, parser_context)
        self.n_scf_iterations = 0
        self.latest_dft_section = None
        self.method_index = None
        self.system_index = None
        self.save_method = False

        # Cache for storing current method settings
        self.method_cache = CacheService(self.parser_context)
        self.method_cache.add("electronic_structure_method", single=False, update=True)

        # Cache for storing current sampling method settings
        self.sampling_method_cache = CacheService(self.parser_context)
        self.sampling_method_cache.add("sampling_method", single=False, update=True)
        self.sampling_method_cache.add("ensemble_type", single=False, update=True)

        # Cache for storing frame sequence information
        self.frame_sequence_cache = CacheService(self.parser_context)
        self.frame_sequence_cache.add("number_of_frames_in_sequence", 0, single=False, update=True)
        self.frame_sequence_cache.add("frame_sequence_local_frames_ref", [], single=False, update=True)
        self.frame_sequence_cache.add("frame_sequence_potential_energy", [], single=False, update=True)
        self.frame_sequence_cache.add("frame_sequence_kinetic_energy", [], single=False, update=True)
        self.frame_sequence_cache.add("frame_sequence_temperature", [], single=False, update=True)
        self.frame_sequence_cache.add("frame_sequence_time", [], single=False, update=True)
        self.frame_sequence_cache.add("frame_sequence_to_sampling_ref", single=False, update=True)

        # Cache for storing system information
        self.system_cache = CacheService(self.parser_context)
        self.system_cache.add("current_positions", single=False, update=False)
        self.system_cache.add("current_labels", single=False, update=False)

        #=======================================================================
        # Cache levels
        self.caching_levels.update({
            'x_nwchem_section_geo_opt_module': CachingLevel.Cache,
            'x_nwchem_section_geo_opt_step': CachingLevel.Cache,
            'x_nwchem_section_xc_functional': CachingLevel.Cache,
            'x_nwchem_section_qmd_module': CachingLevel.ForwardAndCache,
            'x_nwchem_section_qmd_step': CachingLevel.ForwardAndCache,
        })

        #=======================================================================
        # Main Structure
        self.root_matcher = SM("",
            forwardMatch=True,
            sections=['section_run'],
            subMatchers=[
                self.input(),
                self.header(),
                self.system(),

                # This repeating submatcher supports multiple different tasks
                # within one run
                SM("(\s+NWChem DFT Module)|(\s+NWChem Geometry Optimization)|(\s+NWChem QMD Module)",
                    repeats=True,
                    forwardMatch=True,
                    subFlags=SM.SubFlags.Unordered,
                    subMatchers=[
                        self.energy_force_task(),
                        self.geo_opt_module(),
                        self.dft_gaussian_md_task(),
                    ]
                ),
            ]
        )

    def input(self):
        """Returns the simplematcher that parses the NWChem input
        """
        return SM( re.escape("============================== echo of input deck =============================="),
            endReStr=re.escape("================================================================================"),
        )

    def system(self):
        return SM( "                                NWChem Input Module",
            subMatchers=[
                self.geometry(),
                SM( r"  No\.       Tag          Charge          X              Y              Z"),
                SM( re.escape(r" ---- ---------------- ---------- -------------- -------------- --------------"),
                    adHoc=self.adHoc_atoms,
                ),
            ]
        )

    def energy_force_task(self):
        return SM( "                                 NWChem DFT Module",
            forwardMatch=True,
            sections=["section_single_configuration_calculation", "section_system", "section_method"],
            subMatchers=[
                self.dft_module(dft_on_close=self.save_dft_data, scf_on_close=self.save_scf_data, force_on_close=self.save_force_data),
            ],
        )

    def geo_opt_module(self):
        return SM( "                           NWChem Geometry Optimization",
            sections=["section_method", "section_frame_sequence", "section_sampling_method", "x_nwchem_section_geo_opt_module"],
            # onClose={
                # "section_sampling_method": self.save_geo_opt_sampling_id,
            # },
            subFlags=SM.SubFlags.Sequenced,
            subMatchers=[
                SM( r" maximum gradient threshold         \(gmax\) =\s+(?P<geometry_optimization_threshold_force__forceAu>{})".format(self.regexs.float)),
                SM( r" rms gradient threshold             \(grms\) =\s+{}".format(self.regexs.float)),
                SM( r" maximum cartesian step threshold   \(xmax\) =\s+(?P<geometry_optimization_geometry_change__bohr>{})".format(self.regexs.float)),
                SM( r" rms cartesian step threshold       \(xrms\) =\s+{}".format(self.regexs.float)),
                SM( r" fixed trust radius                \(trust\) =\s+{}".format(self.regexs.float)),
                SM( r" maximum step size to saddle      \(sadstp\) =\s+{}".format(self.regexs.float)),
                SM( r" energy precision                  \(eprec\) =\s+(?P<geometry_optimization_energy_change__hartree>{})".format(self.regexs.float)),
                SM( r" maximum number of steps          \(nptopt\) =\s+{}".format(self.regexs.int)),
                SM( r" initial hessian option           \(inhess\) =\s+{}".format(self.regexs.int)),
                SM( r" line search option               \(linopt\) =\s+{}".format(self.regexs.int)),
                SM( r" hessian update option            \(modupd\) =\s+{}".format(self.regexs.int)),
                SM( r" saddle point option              \(modsad\) =\s+{}".format(self.regexs.int)),
                SM( r" initial eigen-mode to follow     \(moddir\) =\s+{}".format(self.regexs.int)),
                SM( r" initial variable to follow       \(vardir\) =\s+{}".format(self.regexs.int)),
                SM( r" follow first negative mode     \(firstneg\) =\s+{}".format(self.regexs.word)),
                SM( r" apply conjugacy                    \(opcg\) =\s+{}".format(self.regexs.word)),
                SM( r" source of zmatrix                         =\s+{}".format(self.regexs.word)),

                SM("          Step\s+\d+$",
                    endReStr="      Optimization converged",
                    forwardMatch=True,
                    subMatchers=[
                        SM("          Step\s+\d+$",
                            repeats=True,
                            weak=True,
                            forwardMatch=True,
                            sections=["section_single_configuration_calculation", "section_system"],
                            onClose={
                                "section_single_configuration_calculation": self.add_frame_reference
                            },
                            subMatchers=[
                                SM("          Step\s+\d+$",
                                    sections=["x_nwchem_section_geo_opt_step"],
                                    subMatchers=[
                                        self.geometry(),
                                        self.dft_module(dft_on_close=self.save_dft_data, scf_on_close=self.save_scf_data, force_on_close=self.save_force_data),
                                        SM( "[.@] Step       Energy      Delta E   Gmax     Grms     Xrms     Xmax   Walltime"),
                                        SM( "[.@] ---- ---------------- -------- -------- -------- -------- -------- --------"),
                                        SM( "@\s+{0}\s+(?P<x_nwchem_geo_opt_step_energy__hartree>{1})\s+{1}\s+{1}\s+{1}\s+{1}\s+{1}\s+{1}".format(self.regexs.int, self.regexs.float)),
                                        self.dft_module(),
                                    ],
                                ),
                            ]
                        )
                    ]
                )
            ]
        )

    def dft_gaussian_md_task(self):
        return SM( "                                 NWChem QMD Module",
            sections=["section_method", "section_frame_sequence", "section_sampling_method", "x_nwchem_section_qmd_module"],
            subMatchers=[
                SM("                                QMD Run Parameters",
                    sections=["x_nwchem_section_qmd_run_parameters"],
                    subMatchers=[
                        SM("    No. of nuclear steps:\s+(?P<x_nwchem_qmd_number_of_nuclear_steps>{})".format(self.regexs.int)),
                        SM("       Nuclear time step:\s+(?P<x_nwchem_qmd_nuclear_time_step>{})".format(self.regexs.float)),
                        SM("        Target temp\. \(K\):\s+(?P<x_nwchem_qmd_target_temperature>{})".format(self.regexs.float)),
                        SM("              Thermostat:\s+(?P<x_nwchem_qmd_thermostat>{})".format(self.regexs.eol)),
                        SM("                     Tau:\s+(?P<x_nwchem_qmd_tau>{})".format(self.regexs.float)),
                        SM("             Random seed:\s+(?P<x_nwchem_qmd_random_seed>{})".format(self.regexs.int)),
                        SM("      Nuclear integrator:\s+(?P<x_nwchem_qmd_nuclear_integrator>{})".format(self.regexs.eol)),
                        SM("       Current temp. \(K\):\s+(?P<x_nwchem_qmd_initial_temperature__K>{})".format(self.regexs.float)),
                    ]
                ),
                SM("                                 NWChem DFT Module",
                    subMatchers=[
                        SM("                         DFT ENERGY GRADIENTS"),
                    ]
                ),
                SM("                                 NWChem DFT Module",
                    repeats=True,
                    sections=["x_nwchem_section_qmd_step", "section_single_configuration_calculation", "section_system"],
                    onClose={
                        "section_single_configuration_calculation": self.add_frame_reference
                    },
                    subMatchers=[
                        SM("                         DFT ENERGY GRADIENTS"),
                        SM("            QMD Run Information",
                            subMatchers=[
                                SM("  Time elapsed \(fs\) :\s+(?P<x_nwchem_qmd_step_time__fs>{})".format(self.regexs.float)),
                                SM("  Kin. energy \(a\.u\.\):\s+{}\s+(?P<x_nwchem_qmd_step_kinetic_energy__hartree>{})".format(self.regexs.int, self.regexs.float)),
                                SM("  Pot. energy \(a\.u\.\):\s+{}\s+(?P<x_nwchem_qmd_step_potential_energy__hartree>{})".format(self.regexs.int, self.regexs.float)),
                                SM("  Tot. energy \(a\.u\.\):\s+{}\s+(?P<x_nwchem_qmd_step_total_energy__hartree>{})".format(self.regexs.int, self.regexs.float)),
                                SM("  Target temp\. \(K\)  :\s+{}\s+(?P<x_nwchem_qmd_step_target_temperature__K>{})".format(self.regexs.int, self.regexs.float)),
                                SM("  Current temp\. \(K\) :\s+{}\s+(?P<x_nwchem_qmd_step_temperature__K>{})".format(self.regexs.int, self.regexs.float)),
                                SM("  Dipole \(a\.u\.\)     :\s+{0}\s+({1}\s+{1}\s+{1})".format(self.regexs.int, self.regexs.float), startReTransform=self.dipole_transform)
                            ]
                        )
                    ]
                )
            ]
        )

    def geometry(self):
        return SM( r"                         Geometry \"geometry\" -> \"geometry\"",
            sections=["x_nwchem_section_geometry"],
            subMatchers=[
                SM(r"                         ---------------------------------"),
                SM(r" Output coordinates in angstroms \(scale by\s+{}to convert to a\.u\.\)"),
                SM(r"  No\.       Tag          Charge          X              Y              Z"),
                SM(r" ---- ---------------- ---------- -------------- -------------- --------------",
                    adHoc=self.adHoc_atoms),
            ]
        )

    def header(self):
        """Returns the simplematcher that parser the NWChem header
        """
        return SM( "              Northwest Computational Chemistry Package \(NWChem\) (?P<program_version>{})".format(self.regexs.float),
            sections=["x_nwchem_section_start_information"],
            subMatchers=[
                SM( r"\s+hostname\s+= (?P<x_nwchem_run_host_name>{})".format(self.regexs.eol)),
                SM( r"\s+program\s+= (?P<x_nwchem_program_name>{})".format(self.regexs.eol)),
                SM( r"\s+date\s+= (?P<x_nwchem_start_datetime>{})".format(self.regexs.eol)),
                SM( r"\s+compiled\s+= (?P<x_nwchem_compilation_datetime>{})".format(self.regexs.eol)),
                SM( r"\s+compiled\s+= (?P<x_nwchem_compilation_datetime>{})".format(self.regexs.eol)),
                SM( r"\s+source\s+= (?P<x_nwchem_source>{})".format(self.regexs.eol)),
                SM( r"\s+nwchem branch\s+= (?P<x_nwchem_branch>{})".format(self.regexs.eol)),
                SM( r"\s+nwchem revision\s+= (?P<x_nwchem_revision>{})".format(self.regexs.eol)),
                SM( r"\s+ga revision\s+= (?P<x_nwchem_ga_revision>{})".format(self.regexs.eol)),
                SM( r"\s+input\s+= (?P<x_nwchem_input_filename>{})".format(self.regexs.eol)),
                SM( r"\s+prefix\s+= (?P<x_nwchem_input_prefix>{})".format(self.regexs.eol)),
                SM( r"\s+data base\s+= (?P<x_nwchem_db_filename>{})".format(self.regexs.eol)),
                SM( r"\s+status\s+= (?P<x_nwchem_status>{})".format(self.regexs.eol)),
                SM( r"\s+nproc\s+= (?P<x_nwchem_nproc>{})".format(self.regexs.eol)),
                SM( r"\s+time left\s+= (?P<x_nwchem_time_left>{})".format(self.regexs.eol)),
            ]
        )

    def dft_module(self, dft_on_close=None, scf_on_close=None, force_on_close=None):
        return SM( "                                 NWChem DFT Module",
            sections=["x_nwchem_section_dft_module"],
            onClose={"x_nwchem_section_dft_module": dft_on_close},
            subMatchers=[
                SM( r"          No. of atoms     :\s+(?P<x_nwchem_dft_number_of_atoms>{})".format(self.regexs.int)),
                SM( r"          Charge           :\s+(?P<x_nwchem_dft_total_charge>{})".format(self.regexs.int)),
                SM( r"          Spin multiplicity:\s+(?P<x_nwchem_dft_spin_multiplicity>{})".format(self.regexs.int)),
                SM( r"          Maximum number of iterations:\s+(?P<x_nwchem_dft_max_iteration>{})".format(self.regexs.int)),
                SM( r"          Convergence on energy requested:\s+(?P<x_nwchem_dft_scf_threshold_energy_change__hartree>{})".format(self.regexs.float)),
                SM( r"          Convergence on density requested:\s+{}".format(self.regexs.float)),
                SM( r"          Convergence on gradient requested:\s+{}".format(self.regexs.float)),
                SM( r"              XC Information",
                    subFlags=SM.SubFlags.Unordered,
                    subMatchers=[
                        SM("\s+(?P<x_nwchem_xc_functional_shortcut>B3LYP Method XC Potential)"),
                        SM("\s+(?P<x_nwchem_xc_functional_shortcut>PBE0 Method XC Functional)"),
                        SM("\s+(?P<x_nwchem_xc_functional_shortcut>Becke half-and-half Method XC Potential)"),
                        SM("\s+(?P<x_nwchem_xc_functional_shortcut>HCTH120  Method XC Functional)"),
                        SM("\s+(?P<x_nwchem_xc_functional_shortcut>HCTH147  Method XC Functional)"),
                        SM("\s+(?P<x_nwchem_xc_functional_shortcut>HCTH407 Method XC Functional)"),
                        SM("\s+(?P<x_nwchem_xc_functional_name>PerdewBurkeErnzerhof Exchange Functional)\s+(?P<x_nwchem_xc_functional_weight>{})".format(self.regexs.float), sections=["x_nwchem_section_xc_functional"]),
                        SM("\s+(?P<x_nwchem_xc_functional_name>Becke 1988 Exchange Functional)\s+(?P<x_nwchem_xc_functional_weight>{})\s+(?P<x_nwchem_xc_functional_type>{})".format(self.regexs.float, self.regexs.eol), sections=["x_nwchem_section_xc_functional"]),
                        SM("\s+(?P<x_nwchem_xc_functional_name>Lee-Yang-Parr Correlation Functional)\s+(?P<x_nwchem_xc_functional_weight>{})\s+(?P<x_nwchem_xc_functional_type>{})".format(self.regexs.float, self.regexs.eol), sections=["x_nwchem_section_xc_functional"]),
                        SM("\s+(?P<x_nwchem_xc_functional_name>Perdew 1991   Exchange Functional)\s+(?P<x_nwchem_xc_functional_weight>{})\s+(?P<x_nwchem_xc_functional_type>{})".format(self.regexs.float, self.regexs.eol), sections=["x_nwchem_section_xc_functional"]),
                        SM("\s+(?P<x_nwchem_xc_functional_name>Perdew 1991 Correlation Functional)\s+(?P<x_nwchem_xc_functional_weight>{})\s+(?P<x_nwchem_xc_functional_type>{})".format(self.regexs.float, self.regexs.eol), sections=["x_nwchem_section_xc_functional"]),
                        SM("\s+(?P<x_nwchem_xc_functional_name>Perdew 1991 LDA Correlation Functional)\s+(?P<x_nwchem_xc_functional_weight>{})\s+(?P<x_nwchem_xc_functional_type>{})".format(self.regexs.float, self.regexs.eol), sections=["x_nwchem_section_xc_functional"]),
                        SM("\s+(?P<x_nwchem_xc_functional_name>PerdewBurkeErnz. Correlation Functional)\s+(?P<x_nwchem_xc_functional_weight>{})\s+(?P<x_nwchem_xc_functional_type>{})".format(self.regexs.float, self.regexs.eol), sections=["x_nwchem_section_xc_functional"]),
                        SM("\s+(?P<x_nwchem_xc_functional_name>Perdew 1981 Correlation Functional)\s+(?P<x_nwchem_xc_functional_weight>{})\s+(?P<x_nwchem_xc_functional_type>{})".format(self.regexs.float, self.regexs.eol), sections=["x_nwchem_section_xc_functional"]),
                        SM("\s+(?P<x_nwchem_xc_functional_name>Perdew 1986 Correlation Functional)\s+(?P<x_nwchem_xc_functional_weight>{})\s+(?P<x_nwchem_xc_functional_type>{})".format(self.regexs.float, self.regexs.eol), sections=["x_nwchem_section_xc_functional"]),
                        SM("\s+(?P<x_nwchem_xc_functional_name>Perdew 1991 Correlation Functional)\s+(?P<x_nwchem_xc_functional_weight>{})\s+(?P<x_nwchem_xc_functional_type>{})".format(self.regexs.float, self.regexs.eol), sections=["x_nwchem_section_xc_functional"]),
                        SM("\s+(?P<x_nwchem_xc_functional_name>Hartree-Fock \(Exact\) Exchange)\s+(?P<x_nwchem_xc_functional_weight>{})\s+(?P<x_nwchem_xc_functional_type>{})".format(self.regexs.float, self.regexs.eol), sections=["x_nwchem_section_xc_functional"]),
                        SM("\s+(?P<x_nwchem_xc_functional_name>Slater Exchange Functional)\s+(?P<x_nwchem_xc_functional_weight>{})\s+(?P<x_nwchem_xc_functional_type>{})".format(self.regexs.float, self.regexs.eol), sections=["x_nwchem_section_xc_functional"]),
                        SM("\s+(?P<x_nwchem_xc_functional_name>OPTX     Exchange Functional)\s+(?P<x_nwchem_xc_functional_weight>{})\s+(?P<x_nwchem_xc_functional_type>{})".format(self.regexs.float, self.regexs.eol), sections=["x_nwchem_section_xc_functional"]),
                    ],
                ),
                SM( r"   convergence    iter        energy       DeltaE   RMS-Dens  Diis-err    time",
                    sections=["x_nwchem_section_dft_scf"],
                    subMatchers=[
                        SM( r" d=\s+{1},ls={0},diis\s+{1}\s+(?P<x_nwchem_dft_scf_energy__hartree>{0})\s+(?P<x_nwchem_dft_energy_change_scf_iteration__hartree>{0})\s+{0}\s+{0}\s+{0}".format(self.regexs.float, self.regexs.int),
                            sections=["x_nwchem_section_dft_scf_step"],
                            onClose={"x_nwchem_section_dft_scf_step": scf_on_close},
                            repeats=True,
                        )
                    ]
                ),
                SM( r"         Total DFT energy =\s+(?P<x_nwchem_dft_energy_total__hartree>{})".format(self.regexs.float)),
                SM( r"      One electron energy =\s+{}".format(self.regexs.float)),
                SM( r"           Coulomb energy =\s+{}".format(self.regexs.float)),
                SM( r"          Exchange energy =\s+(?P<x_nwchem_dft_energy_X__hartree>{})".format(self.regexs.float)),
                SM( r"       Correlation energy =\s+(?P<x_nwchem_dft_energy_C__hartree>{})".format(self.regexs.float)),
                SM( r" Nuclear repulsion energy =\s+{}".format(self.regexs.float)),
                self.dft_gradient_module(force_on_close),
            ],
        )

    def dft_gradient_module(self, on_close=None):
        return SM( r"                            NWChem DFT Gradient Module",
            sections=["x_nwchem_section_dft_gradient"],
            onClose={"x_nwchem_section_dft_gradient": on_close},
            subMatchers=[
                SM( r"                         DFT ENERGY GRADIENTS"),
                SM( r"    atom               coordinates                        gradient"),
                SM( r"                 x          y          z           x          y          z",
                    adHoc=self.adHoc_forces,
                ),
            ],
        )

    #=======================================================================
    # onClose triggers
    def onClose_section_run(self, backend, gIndex, section):
        backend.addValue("program_name", "NWChem")
        backend.addValue("program_basis_set_type", "gaussians+plane_waves")

    def onClose_section_single_configuration_calculation(self, backend, gIndex, section):
        backend.addValue("single_configuration_to_calculation_method_ref", self.method_index)
        backend.addValue("single_configuration_calculation_to_system_ref", self.system_index)

    def onClose_x_nwchem_section_dft_module(self, backend, gIndex, section):
        self.method_cache["electronic_structure_method"] = "DFT"

    def onClose_section_method(self, backend, gIndex, section):
        self.method_cache.addValue("electronic_structure_method")
        self.method_cache.clear()

    def onClose_section_sampling_method(self, backend, gIndex, section):
        self.sampling_method_cache.addValue("sampling_method")
        if self.sampling_method_cache["ensemble_type"] is not None:
            self.sampling_method_cache.addValue("ensemble_type")
        self.sampling_method_cache.clear()

    def onClose_section_system(self, backend, gIndex, section):
        self.system_cache.addArrayValues("atom_positions", "current_positions", unit="angstrom")
        self.system_cache.addArrayValues("atom_labels", "current_labels")
        self.system_index = gIndex

    def onClose_section_frame_sequence(self, backend, gIndex, section):
        self.frame_sequence_cache.addValue("number_of_frames_in_sequence")
        self.frame_sequence_cache.addArrayValues("frame_sequence_local_frames_ref")
        self.frame_sequence_cache.addValue("frame_sequence_to_sampling_ref")

        potential_energy = np.array(self.frame_sequence_cache["frame_sequence_potential_energy"])
        if potential_energy.size != 0:
            backend.addArrayValues("frame_sequence_potential_energy", potential_energy)
            backend.addArrayValues("frame_sequence_potential_energy_stats", np.array([potential_energy.mean(), potential_energy.std()]))

        kin_energy = np.array(self.frame_sequence_cache["frame_sequence_kinetic_energy"])
        if kin_energy.size != 0:
            backend.addArrayValues("frame_sequence_kinetic_energy", kin_energy)
            backend.addArrayValues("frame_sequence_kinetic_energy_stats", np.array([kin_energy.mean(), kin_energy.std()]))

        temp = np.array(self.frame_sequence_cache["frame_sequence_temperature"])
        if temp.size != 0:
            backend.addArrayValues("frame_sequence_temperature", temp)
            backend.addArrayValues("frame_sequence_temperature_stats", np.array([temp.mean(), temp.std()]))

        time = np.array(self.frame_sequence_cache["frame_sequence_time"])
        if time.size != 0:
            backend.addArrayValues("frame_sequence_time", time)

        self.frame_sequence_cache.clear()

    def onClose_x_nwchem_section_qmd_run_parameters(self, backend, gIndex, section):
        thermostat = section.get_latest_value("x_nwchem_qmd_thermostat")
        ensemble = None
        if thermostat == "svr":
            ensemble = "NVT"
        self.sampling_method_cache["ensemble_type"] = ensemble

    def onClose_x_nwchem_section_qmd_step(self, backend, gIndex, section):
        self.method_cache["electronic_structure_method"] = "DFT"
        self.frame_sequence_cache["number_of_frames_in_sequence"] += 1

        potential_energy = section.get_latest_value("x_nwchem_qmd_step_potential_energy")
        self.frame_sequence_cache["frame_sequence_potential_energy"].append(potential_energy)

        kin_energy = section.get_latest_value("x_nwchem_qmd_step_kinetic_energy")
        self.frame_sequence_cache["frame_sequence_kinetic_energy"].append(kin_energy)

        temp = section.get_latest_value("x_nwchem_qmd_step_temperature")
        self.frame_sequence_cache["frame_sequence_temperature"].append(temp)

        time = section.get_latest_value("x_nwchem_qmd_step_time")
        self.frame_sequence_cache["frame_sequence_time"].append(time)

    def onClose_x_nwchem_section_geo_opt_step(self, backend, gIndex, section):
        self.frame_sequence_cache["number_of_frames_in_sequence"] += 1
        pot_ener = section.get_latest_value("x_nwchem_geo_opt_step_energy")
        self.frame_sequence_cache["frame_sequence_potential_energy"].append(pot_ener)

    #=======================================================================
    # onOpen triggers
    def onOpen_section_method(self, backend, gIndex, section):
        self.method_index = gIndex
        self.save_method = True

    def onOpen_section_sampling_method(self, backend, gIndex, section):
        self.frame_sequence_cache["frame_sequence_to_sampling_ref"] = gIndex

    def onOpen_x_nwchem_section_qmd_module(self, backend, gIndex, section):
        self.sampling_method_cache["sampling_method"] = "molecular_dynamics"

    def onOpen_x_nwchem_section_geo_opt_module(self, backend, gIndex, section):
        self.sampling_method_cache["sampling_method"] = "geometry_optimization"

    def adHoc_forces(self, parser):
        # Define the regex that extracts the information
        regex_string = r"\s+({0})\s+({1})\s+({2})\s+({2})\s+({2})\s+({2})\s+({2})\s+({2})".format(self.regexs.int, self.regexs.word, self.regexs.float)
        regex_compiled = re.compile(regex_string)

        match = True
        forces = []

        while match:
            line = parser.fIn.readline()
            result = regex_compiled.match(line)

            if result:
                match = True
                force = [float(x) for x in result.groups()[5:8]]
                forces.append(force)
            else:
                match = False
        forces = -np.array(forces)

        # If anything found, push the results to the correct section
        if len(forces) != 0:
            self.backend.addArrayValues("x_nwchem_dft_forces", forces, unit="forceAu")

    def adHoc_atoms(self, parser):
        # Define the regex that extracts the information
        regex_string = r"\s+({0})\s+({1})\s+({2})\s+({2})\s+({2})\s+({2})".format(self.regexs.int, self.regexs.word, self.regexs.float)
        regex_compiled = re.compile(regex_string)

        match = True
        coordinates = []
        labels = []

        while match:
            line = parser.fIn.readline()
            result = regex_compiled.match(line)

            if result:
                match = True
                results = result.groups()
                label = results[1]
                labels.append(label)
                coordinate = [float(x) for x in results[3:6]]
                coordinates.append(coordinate)
            else:
                match = False
        coordinates = np.array(coordinates)
        labels = np.array(labels)

        # If anything found, push the results to the correct section
        if len(coordinates) != 0:
            self.system_cache["current_positions"] = coordinates
            self.system_cache["current_labels"] = labels

    #=======================================================================
    # SimpleMatcher specific onClose functions
    def save_dft_data(self, backend, gIndex, section):
        section.add_latest_value("x_nwchem_dft_energy_total", "energy_total")
        section.add_latest_value("x_nwchem_dft_energy_X", "energy_X")
        section.add_latest_value("x_nwchem_dft_energy_C", "energy_C")
        section.add_latest_value("x_nwchem_dft_number_of_atoms", "number_of_atoms")
        backend.addValue("number_of_scf_iterations", self.n_scf_iterations)

        # If a geo opt has just been started, save the general settings
        if self.save_method:

            # Basic settings
            section.add_latest_value("x_nwchem_dft_spin_multiplicity", "spin_target_multiplicity")
            section.add_latest_value("x_nwchem_dft_total_charge", "total_charge")
            section.add_latest_value("x_nwchem_dft_max_iteration", "scf_max_iteration")
            section.add_latest_value("x_nwchem_dft_scf_threshold_energy_change", "scf_threshold_energy_change")

            # XC settings
            class XCFunctional(object):
                def __init__(self, name, weight, locality=None):
                    self.name = name
                    self.weight = weight
                    self.locality = locality
                    self.piece = False

                def __eq__(self, other):
                    if isinstance(other, self.__class__):
                        return self.__dict__ == other.__dict__
                    else:
                        return False

                def get_key(self):
                    return "{}_{}_{}".format(self.name, self.weight, self.locality)

            xc_final_list = []

            # Check if shortcut was defined
            shortcut = section.get_latest_value("x_nwchem_xc_functional_shortcut")
            if shortcut:
                shortcut_map = {
                    "B3LYP Method XC Potential": "HYB_GGA_XC_B3LYP",
                    "PBE0 Method XC Functional": "HYB_GGA_XC_PBEH",
                    "Becke half-and-half Method XC Potential": "HYB_GGA_XC_BHANDH",
                    "HCTH120  Method XC Functional": "GGA_XC_HCTH_120",
                    "HCTH147  Method XC Functional": "GGA_XC_HCTH_147",
                    "HCTH407 Method XC Functional": "GGA_XC_HCTH_407",
                }
                norm_name = shortcut_map.get(shortcut)
                if norm_name:
                    xc_final_list.append(XCFunctional(norm_name, 1.0))
            else:
                # Check if any combination with a more generic name is present
                functionals = section["x_nwchem_section_xc_functional"]
                if functionals:
                    xc_info = {}
                    for functional in functionals:
                        name = functional.get_latest_value("x_nwchem_xc_functional_name")
                        weight = functional.get_latest_value("x_nwchem_xc_functional_weight")
                        locality = functional.get_latest_value("x_nwchem_xc_functional_type")
                        if locality is not None:
                            locality = locality.strip()
                        xc = XCFunctional(name, weight, locality)
                        xc_info[xc.get_key()] = xc

                    combinations = {
                        "GGA_X_OPTX": (
                            XCFunctional("OPTX     Exchange Functional", "1.432", "non-local"),
                            XCFunctional("Slater Exchange Functional", "1.052", "local")
                        ),
                        "GGA_C_PBE": (
                            XCFunctional("Perdew 1991 LDA Correlation Functional", "1.0", "local"),
                            XCFunctional("PerdewBurkeErnz. Correlation Functional", "1.0", "non-local")
                        ),
                        "GGA_C_P86": (
                            XCFunctional("Perdew 1981 Correlation Functional", "1.0", "local"),
                            XCFunctional("Perdew 1986 Correlation Functional", "1.0", "non-local")
                        ),
                        "GGA_C_PW91": (
                            XCFunctional("Perdew 1991 Correlation Functional", "1.0", "non-local"),
                            XCFunctional("Perdew 1991 LDA Correlation Functional", "1.0", "local")
                        ),
                    }
                    for name, parts in combinations.items():
                        combination_found = True
                        for part in parts:
                            if part.get_key() not in xc_info:
                                combination_found = False

                        if combination_found:
                            for part in parts:
                                xc_info[part.get_key()].piece = True
                            xc = XCFunctional(name, 1.0)
                            xc_final_list.append(xc)

                    # Gather the pieces that were not part of any bigger
                    # combination
                    for xc in xc_info.values():
                        if not xc.piece:
                            component_map = {
                                "PerdewBurkeErnzerhof Exchange Functional": "GGA_X_PBE",
                                "Becke 1988 Exchange Functional": "GGA_X_B88",
                                "Lee-Yang-Parr Correlation Functional": "GGA_C_LYP",
                                "Perdew 1986 Correlation Functional": "GGA_C_P86",
                                "Perdew 1991 Correlation Functional": "GGA_C_PW91",
                                "Perdew 1991   Exchange Functional": "GGA_X_PW91",
                                "Hartree-Fock \(Exact\) Exchange": "HF_X",
                            }
                            name = xc.name
                            locality = xc.locality
                            weight = xc.weight
                            norm_name = component_map.get(name)
                            if norm_name and (locality is None or locality == ""):

                                id_xc = backend.openSection("section_XC_functionals")
                                backend.addValue("XC_functional_name", norm_name)
                                if weight is not None:
                                    backend.addValue("XC_functional_weight", weight)
                                backend.closeSection("section_XC_functionals", id_xc)
                                xc = XCFunctional(norm_name, weight)
                                xc_final_list.append(xc)

            # Create the summary string
            xc_final_list.sort(key=lambda x: x.name)
            xc_summary = ""
            for i_xc, xc in enumerate(xc_final_list):
                if i_xc != 0:
                    xc_summary += "+"
                xc_summary += "{}*{}".format(xc.weight, xc.name)
            if xc_summary is not "":
                self.backend.addValue("XC_functional", xc_summary)

            self.save_method = False

        self.n_scf_iterations = 0

    def save_scf_data(self, backend, gIndex, section):
        self.n_scf_iterations += 1
        scf_id = backend.openSection("section_scf_iteration")
        section.add_latest_value("x_nwchem_dft_scf_energy", "energy_total_scf_iteration")
        section.add_latest_value("x_nwchem_dft_energy_change_scf_iteration", "energy_change_scf_iteration")
        backend.closeSection("section_scf_iteration", scf_id)

    def save_force_data(self, backend, gIndex, section):
        section.add_latest_array_values("x_nwchem_dft_forces", "atom_forces")

    def save_geo_opt_sampling_id(self, backend, gIndex, section):
        backend.addValue("frame_sequence_to_sampling_ref", gIndex)

    def add_frame_reference(self, backend, gIndex, section):
        self.frame_sequence_cache["frame_sequence_local_frames_ref"].append(gIndex)

    #=======================================================================
    # Start match transforms
    def dipole_transform(self, backend, groups):
        dipole = groups[0]
        components = np.array([float(x) for x in dipole.split()])
        backend.addArrayValues("x_nwchem_qmd_step_dipole", components)

    #=======================================================================
    # Misc
    def debug_end(self):
        def wrapper():
            print("DEBUG END")
        return wrapper

    def debug_close(self):
        def wrapper(backend, gIndex, section):
            print("DEBUG CLOSE")
        return wrapper
