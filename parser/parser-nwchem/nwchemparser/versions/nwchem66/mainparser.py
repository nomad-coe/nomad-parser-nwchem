from __future__ import absolute_import
from nomadcore.simple_parser import SimpleMatcher as SM
from nomadcore.caching_backend import CachingLevel
from nomadcore.baseclasses import MainHierarchicalParser
import re
import logging
import numpy as np
LOGGER = logging.getLogger("nomad")


#===============================================================================
class NWChemMainParser(MainHierarchicalParser):
    """The main parser class that is called for all run types. Parses the CPMD
    output file.
    """
    def __init__(self, file_path, parser_context):
        """
        """
        super(NWChemMainParser, self).__init__(file_path, parser_context)
        self.n_scf_iterations = 0
        self.latest_dft_section = None
        self.frame_sequence_local_frames_ref = []
        self.method_index = None
        self.system_index = None

        #=======================================================================
        # Cache levels
        self.caching_levels.update({
            'x_nwchem_section_geo_opt_task': CachingLevel.Cache,
            'x_nwchem_section_geo_opt_step': CachingLevel.Cache,
        })

        #=======================================================================
        # Globally cached values
        self.cache_service.add("current_positions", single=False, update=False)
        self.cache_service.add("current_labels", single=False, update=False)

        #=======================================================================
        # Main Structure
        self.root_matcher = SM("",
            forwardMatch=True,
            sections=['section_run'],
            subMatchers=[
                self.header(),
                self.input_module(),

                self.energy_force_task(),
                self.geo_opt_task(),

                # SM("(?:\s+NWChem DFT Module)|(?:\s+NWChem Geometry Optimization)",
                    # repeats=True,
                    # # weak=True,
                    # # forwardMatch=True,
                    # # subFlags=SM.SubFlags.Unordered,
                    # # subMatchers=[
                        # # # self.energy_force_task(),
                        # # self.geo_opt_task(),
                    # # ]
                # ),
            ]
        )

    def input(self):
        """Returns the simplematcher that parses the NWChem input
        """
        return SM( re.escape("============================== echo of input deck =============================="),
            endReStr=re.escape("================================================================================"),
        )

    def input_module(self):
        return SM( "                                NWChem Input Module",
            subMatchers=[
                self.geometry(),
                SM( r"  No\.       Tag          Charge          X              Y              Z"),
                SM( re.escape(r" ---- ---------------- ---------- -------------- -------------- --------------"),
                    adHoc=self.adHoc_atoms(),
                ),
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
            sections=["x_nwchem_section_dft"],
            onClose={"x_nwchem_section_dft": dft_on_close},
            subMatchers=[
                SM( r"          No. of atoms     :\s+(?P<x_nwchem_dft_number_of_atoms>{})".format(self.regexs.int)),
                SM( r"          Charge           :\s+(?P<x_nwchem_dft_total_charge>{})".format(self.regexs.int)),
                SM( r"          Spin multiplicity:\s+(?P<x_nwchem_dft_spin_multiplicity>{})".format(self.regexs.int)),
                SM( r"          Maximum number of iterations:\s+(?P<x_nwchem_dft_max_iteration>{})".format(self.regexs.int)),
                SM( r"          Convergence on energy requested:\s+(?P<x_nwchem_dft_scf_threshold_energy_change__hartree>{})".format(self.regexs.float)),
                SM( r"          Convergence on density requested:\s+{}".format(self.regexs.float)),
                SM( r"          Convergence on gradient requested:\s+{}".format(self.regexs.float)),
                SM( r"              XC Information",
                    adHoc=self.adHoc_xc_functionals()
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
                    adHoc=self.adHoc_forces(),
                ),
            ],
        )

    def geometry(self):
        return SM( r"                         Geometry \"geometry\" -> \"geometry\"",
            sections=["x_nwchem_section_geometry"],
            subMatchers=[
                SM(r"                         ---------------------------------"),
                SM(r" Output coordinates in angstroms \(scale by\s+{}to convert to a\.u\.\)"),
                SM(r"  No\.       Tag          Charge          X              Y              Z"),
                SM(r" ---- ---------------- ---------- -------------- -------------- --------------",
                    adHoc=self.adHoc_atoms()),
            ]
        )

    def energy_force_task(self):
        return SM( "                                 NWChem DFT Module",
            forwardMatch=True,
            sections=["section_single_configuration_calculation", "section_system", "section_method", "x_nwchem_section_dft_energy_force_task"],
            # onClose={
                # "section_single_configuration_calculation": self.close_energy_force_single_configuration_calculation()
            # },
            # onOpen={
                # "section_system": self.open_energy_force_section_system()
            # },
            subMatchers=[
                self.dft_module(dft_on_close=self.save_dft_data(), scf_on_close=self.save_scf_data(), force_on_close=self.save_force_data()),
            ],
        )

    def geo_opt_task(self):
        return SM( "                           NWChem Geometry Optimization",
            sections=["section_method", "section_frame_sequence", "section_sampling_method", "x_nwchem_section_geo_opt_task"],
            onClose={
                "section_sampling_method": self.save_geo_opt_sampling_id(),
                "section_frame_sequence": self.save_local_frames_ref(),
            },
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
                                "section_single_configuration_calculation": self.close_geo_opt_single_configuration_calculation()
                            },
                            subMatchers=[
                                SM("          Step\s+\d+$",
                                    sections=["x_nwchem_section_geo_opt_step"],
                                    subMatchers=[
                                        self.geometry(),
                                        self.dft_module(dft_on_close=self.save_dft_data(), scf_on_close=self.save_scf_data(), force_on_close=self.save_force_data()),
                                        SM( "[.@] Step       Energy      Delta E   Gmax     Grms     Xrms     Xmax   Walltime"),
                                        SM( "[.@] ---- ---------------- -------- -------- -------- -------- -------- --------"),
                                        SM( "@\s+{0}\s+(?P<x_nwchem_geo_opt_step_energy>{1})\s+{1}\s+{1}\s+{1}\s+{1}\s+{1}\s+{1}".format(self.regexs.int, self.regexs.float)),
                                        self.dft_module(),
                                    ],
                                ),
                            ]
                        )
                    ]
                )
            ]
        )

    #=======================================================================
    # onClose triggers
    def onClose_section_run(self, backend, gIndex, section):
        backend.addValue("program_name", "NWChem")

    def onClose_section_single_configuration_calculation(self, backend, gIndex, section):
        backend.addValue("single_configuration_to_calculation_method_ref", self.method_index)
        backend.addValue("single_configuration_calculation_to_system_ref", self.system_index)

    def onClose_x_nwchem_section_dft(self, backend, gIndex, section):
        backend.addValue("electronic_structure_method", "DFT")

    def onClose_x_nwchem_section_geo_opt_task(self, backend, gIndex, section):
        steps = section["x_nwchem_section_geo_opt_step"]
        if steps:
            n_steps = len(steps)
            backend.addValue("number_of_frames_in_sequence", n_steps)
            pot_eners = []
            for step in steps:
                pot_ener = step.get_latest_value("x_nwchem_geo_opt_step_energy")
                if pot_ener:
                    pot_eners.append(pot_ener)

            pot_eners = np.array(pot_eners)
            backend.addArrayValues("frame_sequence_potential_energy", pot_eners, unit="hartree")
            backend.addArrayValues("frame_sequence_potential_energy_stats", np.array([pot_eners.mean(), pot_eners.std()]), unit="hartree")

        # Sampling method
        backend.addValue("sampling_method", "geometry_optimization")

    def onClose_section_system(self, backend, gIndex, section):
        self.cache_service.addArrayValues("atom_positions", "current_positions", unit="angstrom")
        self.cache_service.addArrayValues("atom_labels", "current_labels")
        self.system_index = gIndex

    #=======================================================================
    # onOpen triggers
    def onOpen_section_method(self, backend, gIndex, section):
        self.method_index = gIndex

    #=======================================================================
    # adHoc
    def adHoc_xc_functionals(self):
        def wrapper(parser):
            pass
        return wrapper

    def adHoc_forces(self):
        def wrapper(parser):
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

        return wrapper

    def adHoc_atoms(self):
        def wrapper(parser):
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
                self.cache_service["current_positions"] = coordinates
                self.cache_service["current_labels"] = labels

        return wrapper

    #=======================================================================
    # SimpleMatcher specific onClose functions
    def save_dft_data(self):
        def wrapper(backend, gIndex, section):
            section.add_latest_value("x_nwchem_dft_energy_total", "energy_total")
            section.add_latest_value("x_nwchem_dft_energy_X", "energy_X")
            section.add_latest_value("x_nwchem_dft_energy_C", "energy_C")
            section.add_latest_value("x_nwchem_dft_spin_multiplicity", "spin_target_multiplicity")
            section.add_latest_value("x_nwchem_dft_number_of_atoms", "number_of_atoms")
            section.add_latest_value("x_nwchem_dft_total_charge", "total_charge")
            section.add_latest_value("x_nwchem_dft_max_iteration", "scf_max_iteration")
            section.add_latest_value("x_nwchem_dft_scf_threshold_energy_change", "scf_threshold_energy_change")
            backend.addValue("number_of_scf_iterations", self.n_scf_iterations)
            self.n_scf_iterations = 0
        return wrapper

    def save_scf_data(self):
        def wrapper(backend, gIndex, section):
            self.n_scf_iterations += 1
            scf_id = backend.openSection("section_scf_iteration")
            section.add_latest_value("x_nwchem_dft_scf_energy", "energy_total_scf_iteration")
            section.add_latest_value("x_nwchem_dft_energy_change_scf_iteration", "energy_change_scf_iteration")
            backend.closeSection("section_scf_iteration", scf_id)
        return wrapper

    def save_force_data(self):
        def wrapper(backend, gIndex, section):
            section.add_latest_array_values("x_nwchem_dft_forces", "atom_forces")
        return wrapper

    def save_geo_opt_sampling_id(self):
        def wrapper(backend, gIndex, section):
            backend.addValue("frame_sequence_to_sampling_ref", gIndex)
        return wrapper

    def close_geo_opt_single_configuration_calculation(self):
        def wrapper(backend, gIndex, section):
            self.frame_sequence_local_frames_ref.append(gIndex)
        return wrapper

    def save_local_frames_ref(self):
        def wrapper(backend, gIndex, section):
            backend.addArrayValues("frame_sequence_local_frames_ref", np.array(self.frame_sequence_local_frames_ref))
            self.frame_sequence_local_frames_ref = []
        return wrapper
