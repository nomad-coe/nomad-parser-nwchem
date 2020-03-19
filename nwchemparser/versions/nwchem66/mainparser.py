# Copyright 2016-2018 Lauri Himanen, Fawzi Mohamed
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

from __future__ import absolute_import
from nomadcore.simple_parser import SimpleMatcher as SM
from nomadcore.caching_backend import CachingLevel
from nomadcore.baseclasses import MainHierarchicalParser, CacheService
import re
import logging
import numpy as np
LOGGER = logging.getLogger("nomad")


class NWChemMainParser(MainHierarchicalParser):
    """The main parser class that is called for all run types. Parses the NWChem
    output file.
    """
    def __init__(self, parser_context):
        """
        """
        super(NWChemMainParser, self).__init__(parser_context)

        # Cache for storing current method settings
        self.method_cache = CacheService(self.parser_context)
        self.method_cache.add("single_configuration_to_calculation_method_ref", single=False, update=False)
        self.method_cache.add("program_basis_set_type", single=True, update=False)

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
        self.system_cache.add("initial_positions", single=False, update=False)
        self.system_cache.add("atom_positions", single=False, update=True)
        self.system_cache.add("atom_labels", single=False, update=False)
        self.system_cache.add("configuration_periodic_dimensions", single=False, update=False)

        # Cache for storing scc information
        self.scc_cache = CacheService(self.parser_context)
        self.scc_cache.add("atom_forces", single=False, update=True)
        self.scc_cache.add("number_of_scf_iterations", 0, single=False, update=True)
        self.scc_cache.add("single_configuration_calculation_to_system_ref", single=False, update=True)

        #=======================================================================
        # Cache levels
        self.caching_levels.update({
            'x_nwchem_section_geo_opt_module': CachingLevel.Cache,
            'x_nwchem_section_geo_opt_step': CachingLevel.Cache,
            'x_nwchem_section_xc_functional': CachingLevel.Cache,
            'x_nwchem_section_qmd_module': CachingLevel.ForwardAndCache,
            'x_nwchem_section_qmd_step': CachingLevel.ForwardAndCache,
            'x_nwchem_section_xc_part': CachingLevel.ForwardAndCache,
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
                SM("(?:\s+NWChem DFT Module)|(?:\s+NWChem Geometry Optimization)|(?:\s+NWChem QMD Module)|(?:\s+\*               NWPW PSPW Calculation              \*)",
                    repeats=True,
                    forwardMatch=True,
                    subFlags=SM.SubFlags.Unordered,
                    subMatchers=[
                        self.energy_force_gaussian_task(),
                        self.energy_force_pw_task(),
                        self.gaussian_geo_opt_module(),
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

    def energy_force_gaussian_task(self):
        return SM( "                                 NWChem DFT Module",
            startReAction=self.set_gaussian_basis,
            forwardMatch=True,
            subMatchers=[
                self.dft_calculation_full(),
            ],

        )

    def energy_force_pw_task(self):
        return SM( "          \*               NWPW PSPW Calculation              \*",
            sections=["section_single_configuration_calculation", "section_system", "section_method"],
            startReAction=self.set_plane_wave_basis,
            fixedStartValues={
                "electronic_structure_method": "DFT",
            },
            subMatchers=[
                SM(" number of processors used:\s+{}".format(self.regexs.int)),
                SM(" processor grid           :\s+{0} x\s+{0}".format(self.regexs.int)),
                SM(" parallel mapping         :{}".format(self.regexs.eol)),
                SM(" parallel mapping         :{}".format(self.regexs.eol)),
                SM(" number of threads        :{}".format(self.regexs.eol)),
                SM(" parallel io              :{}".format(self.regexs.eol)),
                SM(" options:",
                    subMatchers=[
                        SM("      boundary conditions  =\s+({})".format(self.regexs.word),
                            startReAction=self.transform_periodicity
                        ),
                        SM("      electron spin        =\s+(?P<x_nwchem_electron_spin_restriction>{})".format(self.regexs.word),
                        ),
                        SM("      exchange-correlation =\s+({})".format(self.regexs.eol),
                            startReAction=self.transform_xc
                        ),
                    ]
                ),
                SM(" total charge:\s+({})".format(self.regexs.float),
                    startReAction=self.transform_total_charge,
                ),

                SM(" supercell:",
                    subMatchers=[
                        SM("      cell_name:\s+{}".format(self.regexs.eol)),
                        SM("      lattice:\s+a1=<",
                            forwardMatch=True,
                            adHoc=self.adHoc_lattice("simulation_cell", "angstrom")
                        ),
                        SM("      lattice:\s+a1=<",
                            forwardMatch=True,
                            adHoc=self.adHoc_lattice("x_nwchem_reciprocal_simulation_cell", "angstrom^-1")
                        ),
                        SM("      lattice:\s+a1=<",
                            adHoc=self.adHoc_lattice_parameters
                        ),
                    ]
                ),
            ],
        )

    def gaussian_geo_opt_module(self):
        return SM( "                           NWChem Geometry Optimization",
            startReAction=self.set_gaussian_basis,
            sections=["section_frame_sequence", "section_sampling_method", "x_nwchem_section_geo_opt_module"],
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
                SM( r"          Energy Minimization"),
                SM( r" Names of Z-matrix variables"),
                SM( "          Step\s+\d+$",
                    sections=["x_nwchem_section_geo_opt_step"],
                    subMatchers=[
                        self.geometry(),
                        self.dft_calculation_full(on_close_scc=self.add_frame_reference),
                        SM( "[.@] Step       Energy      Delta E   Gmax     Grms     Xrms     Xmax   Walltime"),
                        SM( "[.@] ---- ---------------- -------- -------- -------- -------- -------- --------"),
                        SM( "@\s+{0}\s+(?P<x_nwchem_geo_opt_step_energy__hartree>{1})\s+{1}\s+{1}\s+{1}\s+{1}\s+{1}\s+{1}".format(self.regexs.int, self.regexs.float)),
                        self.dft_calculation_no_method(),
                    ]
                ),
                SM("          Step\s+\d+$",
                    repeats=True,
                    forwardMatch=True,
                    sections=["x_nwchem_section_geo_opt_step"],
                    subMatchers=[
                        SM("          Step\s+\d+$"),
                        self.geometry(),
                        self.dft_calculation_no_method(on_close_scc=self.add_frame_reference),
                        SM( "[.@] Step       Energy      Delta E   Gmax     Grms     Xrms     Xmax   Walltime"),
                        SM( "[.@] ---- ---------------- -------- -------- -------- -------- -------- --------"),
                        SM( "@\s+{0}\s+(?P<x_nwchem_geo_opt_step_energy__hartree>{1})\s+{1}\s+{1}\s+{1}\s+{1}\s+{1}\s+{1}".format(self.regexs.int, self.regexs.float)),
                        self.dft_calculation_no_method(),
                    ],
                ),
                SM( "\s+Optimization converged",
                    fixedStartValues={
                        "geometry_optimization_converged": True
                    }
                )
            ]
        )

    def dft_gaussian_md_task(self):
        return SM( "                                 NWChem QMD Module",
            startReAction=self.set_gaussian_basis,
            sections=["section_frame_sequence", "section_sampling_method", "x_nwchem_section_qmd_module"],
            subMatchers=[
                SM("                                QMD Run Parameters",
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
                self.dft_calculation_full(),
                SM("                                 NWChem DFT Module",
                    repeats=True,
                    sections=[
                        "x_nwchem_section_qmd_step",
                        "section_single_configuration_calculation",
                        "section_system"
                    ],
                    onClose={
                        "section_single_configuration_calculation": self.add_frame_reference,
                        "section_system": self.push_no_periodicity
                    },
                    subMatchers=[
                        SM("  Caching 1-el integrals"),
                        SM("   Time after variat\. SCF:\s+{}".format(self.regexs.float)),
                        SM("   Time prior to 1st pass:\s+{}".format(self.regexs.float)),
                        SM("         Total DFT energy =\s+(?P<energy_total__hartree>{})".format(self.regexs.float)),
                        SM("      One electron energy =\s+(?P<x_nwchem_energy_one_electron__hartree>{})".format(self.regexs.float)),
                        SM("           Coulomb energy =\s+(?P<x_nwchem_energy_coulomb__hartree>{})".format(self.regexs.float)),
                        SM("    Exchange-Corr\. energy =\s+(?P<energy_XC__hartree>{})".format(self.regexs.float)),
                        SM(" Nuclear repulsion energy =\s+(?P<x_nwchem_energy_nuclear_repulsion__hartree>{})".format(self.regexs.float)),
                        SM(" Numeric\. integr\. density =\s+{}".format(self.regexs.float)),
                        SM("     Total iterative time =\s+(?P<time_calculation__s>{})".format(self.regexs.float)),
                        SM("                         DFT ENERGY GRADIENTS",
                            subMatchers=[
                                SM("    atom               coordinates                        gradient"),
                                SM("                 x          y          z           x          y          z",
                                    adHoc=self.adHoc_forces(save_positions=True)
                                ),
                            ]
                        ),
                        SM("            QMD Run Information",
                            subMatchers=[
                                SM("  Time elapsed \(fs\) :\s+(?P<x_nwchem_qmd_step_time__fs>{})".format(self.regexs.float)),
                                SM("  Kin. energy \(a\.u\.\):\s+{}\s+(?P<x_nwchem_qmd_step_kinetic_energy__hartree>{})".format(self.regexs.int, self.regexs.float)),
                                SM("  Pot. energy \(a\.u\.\):\s+{}\s+(?P<x_nwchem_qmd_step_potential_energy__hartree>{})".format(self.regexs.int, self.regexs.float)),
                                SM("  Tot. energy \(a\.u\.\):\s+{}\s+(?P<x_nwchem_qmd_step_total_energy__hartree>{})".format(self.regexs.int, self.regexs.float)),
                                SM("  Target temp\. \(K\)  :\s+{}\s+(?P<x_nwchem_qmd_step_target_temperature__K>{})".format(self.regexs.int, self.regexs.float)),
                                SM("  Current temp\. \(K\) :\s+{}\s+(?P<x_nwchem_qmd_step_temperature__K>{})".format(self.regexs.int, self.regexs.float)),
                                SM("  Dipole \(a\.u\.\)     :\s+{0}\s+({1}\s+{1}\s+{1})".format(self.regexs.int, self.regexs.float), startReAction=self.transform_dipole)
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
        return SM( "\s+Northwest Computational Chemistry Package \(NWChem\) (?P<program_version>(\d+\.\d+(?:\.\d+)?))",
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

    def dft_calculation_no_method(self, on_close_scc=None):
        """SimpleMatcher that ignores DFT method settings but parses actual
        calculation. Used in MD and geo opt, when the settings have not changed
        to avoid storing redundant data.
        """
        return SM( "                                 NWChem DFT Module",
            sections=[
                "section_single_configuration_calculation",
                "section_system",
            ],
            onClose={
                "section_single_configuration_calculation": on_close_scc,
                "section_system": self.push_no_periodicity
            },
            subMatchers=self.dft_calculation_data(),
        )

    def dft_calculation_full(self, on_close_scc=None):
        subMatchers = [
            SM( r"          No. of atoms     :\s+(?P<number_of_atoms>{})".format(self.regexs.int)),
            SM( r"          Charge           :\s+(?P<total_charge>{})".format(self.regexs.int)),
            SM( r"          Spin multiplicity:\s+(?P<spin_target_multiplicity>{})".format(self.regexs.int)),
            SM( r"          Maximum number of iterations:\s+(?P<scf_max_iteration>{})".format(self.regexs.int)),
            SM( r"          Convergence on energy requested:\s+(?P<scf_threshold_energy_change__hartree>{})".format(self.regexs.float)),
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
                    SM("\s+(?P<x_nwchem_xc_functional_name>PerdewBurkeErnzerhof Exchange Functional)\s+(?P<x_nwchem_xc_functional_weight>{})".format(self.regexs.float), sections=["x_nwchem_section_xc_part"]),
                    SM("\s+(?P<x_nwchem_xc_functional_name>Becke 1988 Exchange Functional)\s+(?P<x_nwchem_xc_functional_weight>{})\s+(?P<x_nwchem_xc_functional_type>{})".format(self.regexs.float, self.regexs.eol), sections=["x_nwchem_section_xc_part"]),
                    SM("\s+(?P<x_nwchem_xc_functional_name>Lee-Yang-Parr Correlation Functional)\s+(?P<x_nwchem_xc_functional_weight>{})\s+(?P<x_nwchem_xc_functional_type>{})".format(self.regexs.float, self.regexs.eol), sections=["x_nwchem_section_xc_part"]),
                    SM("\s+(?P<x_nwchem_xc_functional_name>Perdew 1991   Exchange Functional)\s+(?P<x_nwchem_xc_functional_weight>{})\s+(?P<x_nwchem_xc_functional_type>{})".format(self.regexs.float, self.regexs.eol), sections=["x_nwchem_section_xc_part"]),
                    SM("\s+(?P<x_nwchem_xc_functional_name>Perdew 1991 Correlation Functional)\s+(?P<x_nwchem_xc_functional_weight>{})\s+(?P<x_nwchem_xc_functional_type>{})".format(self.regexs.float, self.regexs.eol), sections=["x_nwchem_section_xc_part"]),
                    SM("\s+(?P<x_nwchem_xc_functional_name>Perdew 1991 LDA Correlation Functional)\s+(?P<x_nwchem_xc_functional_weight>{})\s+(?P<x_nwchem_xc_functional_type>{})".format(self.regexs.float, self.regexs.eol), sections=["x_nwchem_section_xc_part"]),
                    SM("\s+(?P<x_nwchem_xc_functional_name>PerdewBurkeErnz. Correlation Functional)\s+(?P<x_nwchem_xc_functional_weight>{})\s+(?P<x_nwchem_xc_functional_type>{})".format(self.regexs.float, self.regexs.eol), sections=["x_nwchem_section_xc_part"]),
                    SM("\s+(?P<x_nwchem_xc_functional_name>Perdew 1981 Correlation Functional)\s+(?P<x_nwchem_xc_functional_weight>{})\s+(?P<x_nwchem_xc_functional_type>{})".format(self.regexs.float, self.regexs.eol), sections=["x_nwchem_section_xc_part"]),
                    SM("\s+(?P<x_nwchem_xc_functional_name>Perdew 1986 Correlation Functional)\s+(?P<x_nwchem_xc_functional_weight>{})\s+(?P<x_nwchem_xc_functional_type>{})".format(self.regexs.float, self.regexs.eol), sections=["x_nwchem_section_xc_part"]),
                    SM("\s+(?P<x_nwchem_xc_functional_name>Perdew 1991 Correlation Functional)\s+(?P<x_nwchem_xc_functional_weight>{})\s+(?P<x_nwchem_xc_functional_type>{})".format(self.regexs.float, self.regexs.eol), sections=["x_nwchem_section_xc_part"]),
                    SM("\s+(?P<x_nwchem_xc_functional_name>Hartree-Fock \(Exact\) Exchange)\s+(?P<x_nwchem_xc_functional_weight>{})\s+(?P<x_nwchem_xc_functional_type>{})".format(self.regexs.float, self.regexs.eol), sections=["x_nwchem_section_xc_part"]),
                    SM("\s+(?P<x_nwchem_xc_functional_name>Slater Exchange Functional)\s+(?P<x_nwchem_xc_functional_weight>{})\s+(?P<x_nwchem_xc_functional_type>{})".format(self.regexs.float, self.regexs.eol), sections=["x_nwchem_section_xc_part"]),
                    SM("\s+(?P<x_nwchem_xc_functional_name>OPTX     Exchange Functional)\s+(?P<x_nwchem_xc_functional_weight>{})\s+(?P<x_nwchem_xc_functional_type>{})".format(self.regexs.float, self.regexs.eol), sections=["x_nwchem_section_xc_part"]),
                    SM("\s+(?P<x_nwchem_xc_functional_name>TPSS metaGGA Exchange Functional)\s+(?P<x_nwchem_xc_functional_weight>{})".format(self.regexs.float, self.regexs.eol), sections=["x_nwchem_section_xc_part"]),
                    SM("\s+(?P<x_nwchem_xc_functional_name>TPSS03 metaGGA Correlation Functional)\s+(?P<x_nwchem_xc_functional_weight>{})".format(self.regexs.float, self.regexs.eol), sections=["x_nwchem_section_xc_part"]),
                    SM("\s+(?P<x_nwchem_xc_functional_name>VWN V Correlation Functional)\s+(?P<x_nwchem_xc_functional_weight>{})".format(self.regexs.float, self.regexs.eol), sections=["x_nwchem_section_xc_part"]),
                ],
            ),
        ]
        subMatchers.extend(self.dft_calculation_data())

        return SM( "                                 NWChem DFT Module",
            sections=[
                "section_single_configuration_calculation",
                "section_system",
                "section_method"
            ],
            fixedStartValues={"electronic_structure_method": "DFT"},
            onClose={
                "section_single_configuration_calculation": on_close_scc,
                "section_system": self.push_no_periodicity
            },
            subMatchers=subMatchers,
        )

    def dft_calculation_data(self):
        return [
            SM( r"   convergence    iter        energy       DeltaE   RMS-Dens  Diis-err    time",
                subMatchers=[
                    SM( r" d=\s+{1},ls={0},diis\s+{1}\s+(?P<energy_total_scf_iteration__hartree>{0})\s+(?P<energy_change_scf_iteration__hartree>{0})\s+{0}\s+{0}\s+{0}".format(self.regexs.float, self.regexs.int),
                        sections=["section_scf_iteration"],
                        repeats=True,
                    )
                ]
            ),
            SM( r"         Total DFT energy =\s+(?P<energy_total__hartree>{})".format(self.regexs.float)),
            SM( r"      One electron energy =\s+(?P<x_nwchem_energy_one_electron__hartree>{})".format(self.regexs.float)),
            SM( r"           Coulomb energy =\s+(?P<x_nwchem_energy_coulomb__hartree>{})".format(self.regexs.float)),
            SM( r"          Exchange energy =\s+(?P<energy_X__hartree>{})".format(self.regexs.float)),
            SM( r"       Correlation energy =\s+(?P<energy_C__hartree>{})".format(self.regexs.float)),
            SM( r" Nuclear repulsion energy =\s+(?P<x_nwchem_energy_nuclear_repulsion__hartree>{})".format(self.regexs.float)),
            self.dft_gradient_module(),
        ]

    def dft_gradient_module(self):
        return SM( r"                            NWChem DFT Gradient Module",
            subMatchers=[
                SM( r"                         DFT ENERGY GRADIENTS"),
                SM( r"    atom               coordinates                        gradient"),
                SM( r"                 x          y          z           x          y          z",
                    adHoc=self.adHoc_forces(),
                ),
            ],
        )

    #=======================================================================
    # onClose triggers
    def onClose_section_run(self, backend, gIndex, section):
        backend.addValue("program_name", "NWChem")
        self.method_cache.addValue("program_basis_set_type")
        # backend.addValue("program_basis_set_type", "gaussians+plane_waves")

    def onClose_section_single_configuration_calculation(self, backend, gIndex, section):
        self.scc_cache.addValue("single_configuration_calculation_to_system_ref")
        self.method_cache.addValue("single_configuration_to_calculation_method_ref")

        if self.scc_cache["atom_forces"] is not None:
            self.scc_cache.addArrayValues("atom_forces", unit="hartree / bohr")

        if self.scc_cache["number_of_scf_iterations"]:
            self.scc_cache.addValue("number_of_scf_iterations")

        self.scc_cache.clear()

    def onClose_section_method(self, backend, gIndex, section):

        # XC settings
        class XCFunctional(object):
            def __init__(self, name, weight, locality=None):
                self.name = name
                self.weight = weight
                self.locality = locality

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
            functionals = section["x_nwchem_section_xc_part"]
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

                # Check if a predetermined combination of
                # x_nwchem_section_xc_parts are found.
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
                            del xc_info[part.get_key()]
                        xc = XCFunctional(name, 1.0)
                        xc_final_list.append(xc)

                # Gather the pieces that were not part of any bigger
                # combination
                for xc in xc_info.values():
                    component_map = {
                        "PerdewBurkeErnzerhof Exchange Functional": "GGA_X_PBE",
                        "Becke 1988 Exchange Functional": "GGA_X_B88",
                        "Lee-Yang-Parr Correlation Functional": "GGA_C_LYP",
                        "Perdew 1986 Correlation Functional": "GGA_C_P86",
                        "Perdew 1991 Correlation Functional": "GGA_C_PW91",
                        "Perdew 1991   Exchange Functional": "GGA_X_PW91",
                        "Hartree-Fock \(Exact\) Exchange": "HF_X",
                        "TPSS metaGGA Exchange Functional": "MGGA_X_TPSS",
                        "TPSS03 metaGGA Correlation Functional": "MGGA_C_TPSS",
                        "Slater Exchange Functional": "LDA_X",
                        "VWN V Correlation Functional": "LDA_C_VWN",
                    }
                    name = xc.name
                    locality = xc.locality
                    weight = xc.weight
                    norm_name = component_map.get(name)
                    if norm_name:
                        xc = XCFunctional(norm_name, weight)
                        xc_final_list.append(xc)

        # Go throught the list of found functionals and push XC sections and
        # gather the summary string.
        xc_final_list.sort(key=lambda x: x.name)
        xc_summary = ""
        for i_xc, xc in enumerate(xc_final_list):
            if i_xc != 0:
                xc_summary += "+"
            xc_summary += "{}*{}".format(xc.weight, xc.name)

            # Push the XC sections
            id_xc = backend.openSection("section_XC_functionals")
            backend.addValue("XC_functional_name", xc.name)
            if xc.weight is not None:
                backend.addValue("XC_functional_weight", xc.weight)
            backend.closeSection("section_XC_functionals", id_xc)

        # Push the summary string
        if xc_summary is not "":
            self.backend.addValue("XC_functional", xc_summary)

    def onClose_section_sampling_method(self, backend, gIndex, section):

        thermostat = section.get_latest_value("x_nwchem_qmd_thermostat")
        ensemble = None
        if thermostat == "svr":
            ensemble = "NVT"
            self.backend.addValue("ensemble_type", ensemble)

        self.sampling_method_cache.addValue("sampling_method")
        self.sampling_method_cache.clear()

    def onClose_section_system(self, backend, gIndex, section):
        if self.system_cache["atom_positions"] is not None:
            self.system_cache.addArrayValues("atom_positions", unit="bohr")
        elif self.system_cache["initial_positions"] is not None:
            self.system_cache.addArrayValues("atom_positions", "initial_positions", unit="angstrom")
        self.system_cache.addArrayValues("atom_labels")
        self.system_cache.addArrayValues("configuration_periodic_dimensions")

    def onClose_section_frame_sequence(self, backend, gIndex, section):
        self.frame_sequence_cache.addValue("number_of_frames_in_sequence")
        frame_sequence = self.frame_sequence_cache["frame_sequence_local_frames_ref"]
        if frame_sequence:
            frame_sequence = np.array(frame_sequence)
            self.backend.addArrayValues("frame_sequence_local_frames_ref", frame_sequence)
        self.frame_sequence_cache.addValue("frame_sequence_to_sampling_ref")

        potential_energy = self.frame_sequence_cache["frame_sequence_potential_energy"]
        if potential_energy:
            potential_energy = np.array(potential_energy)
            backend.addArrayValues("frame_sequence_potential_energy", potential_energy)
            backend.addArrayValues("frame_sequence_potential_energy_stats", np.array([potential_energy.mean(), potential_energy.std()]))

        kin_energy = self.frame_sequence_cache["frame_sequence_kinetic_energy"]
        if kin_energy:
            kin_energy = np.array(kin_energy)
            backend.addArrayValues("frame_sequence_kinetic_energy", kin_energy)
            backend.addArrayValues("frame_sequence_kinetic_energy_stats", np.array([kin_energy.mean(), kin_energy.std()]))

        temp = self.frame_sequence_cache["frame_sequence_temperature"]
        if temp:
            temp = np.array(temp)
            backend.addArrayValues("frame_sequence_temperature", temp)
            backend.addArrayValues("frame_sequence_temperature_stats", np.array([temp.mean(), temp.std()]))

        time = self.frame_sequence_cache["frame_sequence_time"]
        if time:
            time = np.array(time)
            backend.addArrayValues("frame_sequence_time", time)

        self.frame_sequence_cache.clear()

    def onClose_section_scf_iteration(self, backend, gIndex, section):
        self.scc_cache["number_of_scf_iterations"] += 1

    def onClose_x_nwchem_section_qmd_step(self, backend, gIndex, section):
        self.frame_sequence_cache["number_of_frames_in_sequence"] += 1

        potential_energy = section.get_latest_value("x_nwchem_qmd_step_potential_energy")
        if potential_energy is not None:
            self.frame_sequence_cache["frame_sequence_potential_energy"].append(potential_energy)

        kin_energy = section.get_latest_value("x_nwchem_qmd_step_kinetic_energy")
        if kin_energy is not None:
            self.frame_sequence_cache["frame_sequence_kinetic_energy"].append(kin_energy)

        temp = section.get_latest_value("x_nwchem_qmd_step_temperature")
        if temp is not None:
            self.frame_sequence_cache["frame_sequence_temperature"].append(temp)

        time = section.get_latest_value("x_nwchem_qmd_step_time")
        if time is not None:
            self.frame_sequence_cache["frame_sequence_time"].append(time)

    def onClose_x_nwchem_section_geo_opt_step(self, backend, gIndex, section):
        self.frame_sequence_cache["number_of_frames_in_sequence"] += 1
        pot_ener = section.get_latest_value("x_nwchem_geo_opt_step_energy")
        if pot_ener is not None:
            self.frame_sequence_cache["frame_sequence_potential_energy"].append(pot_ener)

    #=======================================================================
    # onOpen triggers
    def onOpen_section_method(self, backend, gIndex, section):
        self.method_cache["single_configuration_to_calculation_method_ref"] = gIndex

    def onOpen_section_system(self, backend, gIndex, section):
        self.scc_cache["single_configuration_calculation_to_system_ref"] = gIndex

    def onOpen_section_sampling_method(self, backend, gIndex, section):
        self.frame_sequence_cache["frame_sequence_to_sampling_ref"] = gIndex

    def onOpen_x_nwchem_section_qmd_module(self, backend, gIndex, section):
        self.sampling_method_cache["sampling_method"] = "molecular_dynamics"

    def onOpen_x_nwchem_section_geo_opt_module(self, backend, gIndex, section):
        self.sampling_method_cache["sampling_method"] = "geometry_optimization"

    #=======================================================================
    # adHoc
    def adHoc_forces(self, save_positions=False):
        def wrapper(parser):
            match = True
            forces = []
            positions = []

            while match:
                line = parser.fIn.readline()
                if line == "" or line.isspace():
                    match = False
                    break
                components = line.split()
                position = np.array([float(x) for x in components[-6:-3]])
                force = np.array([float(x) for x in components[-3:]])
                forces.append(force)
                positions.append(position)

            forces = -np.array(forces)
            positions = np.array(positions)

            # If anything found, push the results to the correct section
            if forces.size != 0:
                self.scc_cache["atom_forces"] = forces
            if save_positions:
                if positions.size != 0:
                    self.system_cache["atom_positions"] = positions
        return wrapper

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
            self.system_cache["initial_positions"] = coordinates
            self.system_cache["atom_labels"] = labels

    def adHoc_lattice(self, metaname, units):
        def wrapper(parser):
            a = [float(x) for x in parser.fIn.readline().split()[2:5]]
            b = [float(x) for x in parser.fIn.readline().split()[1:4]]
            c = [float(x) for x in parser.fIn.readline().split()[1:4]]

            cell = np.array([a, b, c])
            parser.backend.addArrayValues(metaname, cell, unit=units)
        return wrapper

    def adHoc_lattice_parameters(self, parser):
        pass

    #=======================================================================
    # SimpleMatcher specific onClose
    def save_geo_opt_sampling_id(self, backend, gIndex, section):
        backend.addValue("frame_sequence_to_sampling_ref", gIndex)

    def add_frame_reference(self, backend, gIndex, section):
        self.frame_sequence_cache["frame_sequence_local_frames_ref"].append(gIndex)

    def push_no_periodicity(self, backend, gIndex, section):
        self.system_cache["configuration_periodic_dimensions"] = np.array([False, False, False])

    #=======================================================================
    # Start match transforms
    def transform_dipole(self, backend, groups):
        dipole = groups[0]
        components = np.array([float(x) for x in dipole.split()])
        backend.addArrayValues("x_nwchem_qmd_step_dipole", components)

    def transform_periodicity(self, backend, groups):
        periodicity = groups[0]
        if periodicity == "periodic":
            self.system_cache["configuration_periodic_dimensions"] = np.array([True, True, True])

    def transform_xc(self, backend, groups):
        xc = groups[0]
        if xc == "LDA (Vosko et al) parameterization":
            pass

    def transform_total_charge(self, backend, groups):
        charge = groups[0]
        self.backend.addValue("total_charge", int(float(charge)))

    def set_gaussian_basis(self, backend, groups):
        self.method_cache["program_basis_set_type"] = "gaussians"

    def set_plane_wave_basis(self, backend, groups):
        self.method_cache["program_basis_set_type"] = "plane waves"

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
