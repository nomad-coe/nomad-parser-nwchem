"""
This is a module for unit testing the NWChem parser. The unit tests are run with
a custom backend that outputs the results directly into native python object for
easier and faster analysis.

Each property that has an enumerable list of different possible options is
assigned a new test class, that should ideally test through all the options.

The properties that can have non-enumerable values will be tested only for one
specific case inside a test class that is designed for a certain type of run
(MD, optimization, QM/MM, etc.)
"""
import os
import unittest
import logging
import numpy as np
from nwchemparser import NWChemParser
from nomadcore.unit_conversion.unit_conversion import convert_unit


#===============================================================================
def get_results(folder, metainfo_to_keep=None):
    """Get the given result from the calculation in the given folder by using
    the Analyzer in the nomadtoolkit package. Tries to optimize the parsing by
    giving the metainfo_to_keep argument.

    Args:
        folder: The folder relative to the directory of this script where the
            parsed calculation resides.
        metaname: The quantity to extract.
    """
    dirname = os.path.dirname(__file__)
    filename = os.path.join(dirname, folder, "output.out")
    parser = NWChemParser(filename, None, debug=True, log_level=logging.WARNING)
    results = parser.parse()
    return results


#===============================================================================
def get_result(folder, metaname, optimize=True):
    if optimize:
        results = get_results(folder, None)
    else:
        results = get_results(folder)
    result = results[metaname]
    return result


#===============================================================================
class TestDFTGaussianEnergy(unittest.TestCase):
    """Tests that the parser can handle DFT energy calculations.
    """

    @classmethod
    def setUpClass(cls):
        cls.results = get_results("dft_gaussian/energy", "section_run")
        # cls.results.print_summary()

    def test_program_name(self):
        result = self.results["program_name"]
        self.assertEqual(result, "NWChem")

    def test_configuration_periodic_dimensions(self):
        result = self.results["configuration_periodic_dimensions"]
        self.assertTrue(np.array_equal(result, np.array([False, False, False])))

    def test_program_version(self):
        result = self.results["program_version"]
        self.assertEqual(result, "6.6")

    def test_xc_functional(self):
        result = self.results["XC_functional"]
        self.assertEqual(result, "1.0*MGGA_C_TPSS+1.0*MGGA_X_TPSS")

    def test_atom_labels(self):
        atom_labels = self.results["atom_labels"]
        expected_labels = np.array(["O", "H", "H"])
        self.assertTrue(np.array_equal(atom_labels, expected_labels))

    def test_atom_positions(self):
        atom_position = self.results["atom_positions"]
        expected_position = convert_unit(np.array(
            [
                [0.00000000, 0.00000000, -0.11817375],
                [0.76924532, 0.00000000, 0.47269501],
                [-0.76924532, 0.00000000, 0.47269501],
            ]
        ), "angstrom")
        self.assertTrue(np.array_equal(atom_position, expected_position))

    def test_number_of_atoms(self):
        n_atoms = self.results["number_of_atoms"]
        self.assertEqual(n_atoms, 3)

    def test_total_charge(self):
        charge = self.results["total_charge"]
        self.assertEqual(charge, 0)

    def test_energy_total(self):
        result = self.results["energy_total"]
        expected_result = convert_unit(np.array(-76.436222730188), "hartree")
        self.assertTrue(np.array_equal(result, expected_result))

    def test_energy_x(self):
        result = self.results["energy_X"]
        expected_result = convert_unit(np.array(-9.025345841743), "hartree")
        self.assertTrue(np.array_equal(result, expected_result))

    def test_energy_c(self):
        result = self.results["energy_C"]
        expected_result = convert_unit(np.array(-0.328011552453), "hartree")
        self.assertTrue(np.array_equal(result, expected_result))

    def test_energy_total_scf_iteration(self):
        result = self.results["energy_total_scf_iteration"]
        # Test the first and last energies
        expected_result = convert_unit(np.array(
            [
                [-76.3916403957],
                [-76.4362227302],
            ]), "hartree")
        self.assertTrue(np.array_equal(np.array([[result[0]], [result[-1]]]), expected_result))

    def test_energy_change_scf_iteration(self):
        result = self.results["energy_change_scf_iteration"]
        expected_result = convert_unit(np.array(
            [
                [-8.55E+01],
                [-3.82E-07],
            ]), "hartree")
        self.assertTrue(np.array_equal(np.array([[result[0]], [result[-1]]]), expected_result))

    def test_scf_max_iteration(self):
        result = self.results["scf_max_iteration"]
        self.assertEqual(result, 50)

    def test_scf_threshold_energy_change(self):
        result = self.results["scf_threshold_energy_change"]
        self.assertEqual(result, convert_unit(1.00E-06, "hartree"))

    def test_electronic_structure_method(self):
        result = self.results["electronic_structure_method"]
        self.assertEqual(result, "DFT")

    def test_scf_dft_number_of_iterations(self):
        result = self.results["number_of_scf_iterations"]
        self.assertEqual(result, 6)

    def test_spin_target_multiplicity(self):
        multiplicity = self.results["spin_target_multiplicity"]
        self.assertEqual(multiplicity, 1)

    def test_single_configuration_to_calculation_method_ref(self):
        result = self.results["single_configuration_to_calculation_method_ref"]
        self.assertEqual(result, 0)

    def test_single_configuration_calculation_to_system_description_ref(self):
        result = self.results["single_configuration_calculation_to_system_ref"]
        self.assertEqual(result, 0)

    # def test_single_configuration_calculation_converged(self):
        # result = self.results["single_configuration_calculation_converged"]
        # self.assertTrue(result)

    # def test_section_method_atom_kind(self):
        # kind = self.results["section_method_atom_kind"][0]
        # self.assertEqual(kind["method_atom_kind_atom_number"][0], 1)
        # self.assertEqual(kind["method_atom_kind_label"][0], "H")

    # def test_section_method_basis_set(self):
        # kind = self.results["section_method_basis_set"][0]
        # self.assertEqual(kind["method_basis_set_kind"][0], "wavefunction")
        # self.assertTrue(np.array_equal(kind["mapping_section_method_basis_set_cell_associated"][0], 0))

    # def test_number_of_spin_channels(self):
        # result = self.results["number_of_spin_channels"]
        # self.assertEqual(result, 1)

    # def test_simulation_cell(self):
        # cell = self.results["simulation_cell"]
        # n_vectors = cell.shape[0]
        # n_dim = cell.shape[1]
        # self.assertEqual(n_vectors, 3)
        # self.assertEqual(n_dim, 3)
        # expected_cell = convert_unit(np.array([[15.1178, 0, 0], [0, 15.1178, 0], [0, 0, 15.1178]]), "bohr")
        # self.assertTrue(np.array_equal(cell, expected_cell))

    # def test_basis_set_cell_dependent(self):
        # kind = self.results["basis_set_cell_dependent_kind"]
        # name = self.results["basis_set_cell_dependent_name"]
        # cutoff = self.results["basis_set_planewave_cutoff"]

        # self.assertEqual(kind, "plane_waves")
        # self.assertEqual(name, "PW_70.0")
        # self.assertEqual(cutoff, convert_unit(70.00000, "rydberg"))


#===============================================================================
class TestDFTGaussianForce(unittest.TestCase):
    """Tests that the parser can handle DFT force calculations.
    """
    @classmethod
    def setUpClass(cls):
        cls.results = get_results("dft_gaussian/force", "section_run")

    def test_configuration_periodic_dimensions(self):
        result = self.results["configuration_periodic_dimensions"]
        self.assertTrue(np.array_equal(result, np.array([False, False, False])))

    def test_electronic_structure_method(self):
        result = self.results["electronic_structure_method"]
        self.assertEqual(result, "DFT")

    def test_xc_functional(self):
        result = self.results["XC_functional"]
        self.assertEqual(result, "1.0*MGGA_C_TPSS+1.0*MGGA_X_TPSS")

    def test_atom_forces(self):
        result = self.results["atom_forces"]
        expected_result = convert_unit(
            -np.array([
                [0.000000, 0.000000, -0.000037],
                [0.000006, 0.000000, 0.000018],
                [-0.000006, 0.000000, 0.000018],
            ]),
            "hartree/bohr"
        )
        self.assertTrue(np.array_equal(result, expected_result))


#===============================================================================
class TestDFTGaussianGeoOpt(unittest.TestCase):
    """Tests that the parser can handle DFT geometry optimizations.
    """
    @classmethod
    def setUpClass(cls):
        cls.results = get_results("dft_gaussian/geo_opt", "section_run")

    def test_configuration_periodic_dimensions(self):
        result = self.results["configuration_periodic_dimensions"][0]
        self.assertTrue(np.array_equal(result, np.array([False, False, False])))

    def test_xc_functional(self):
        result = self.results["XC_functional"]
        self.assertEqual(result, "1.0*MGGA_C_TPSS+1.0*MGGA_X_TPSS")

    def test_electronic_structure_method(self):
        result = self.results["electronic_structure_method"]
        self.assertEqual(result, "DFT")

    def test_frame_sequence(self):
        sequence = self.results["section_frame_sequence"][0]

        # Number of frames
        n_frames = self.results["number_of_frames_in_sequence"]
        self.assertEqual(n_frames, 4)

        # Potential energy
        pot_ener = sequence["frame_sequence_potential_energy"]
        expected_pot_ener = convert_unit(
            np.array([
                -76.42941861,
                -76.43609119,
                -76.43622176,
                -76.43622273,
            ]),
            "hartree"
        )
        self.assertTrue(np.array_equal(pot_ener, expected_pot_ener))

        # Test positions
        positions = self.results["atom_positions"]
        expected_pos = convert_unit(
            np.array([
                [
                    [0.00000000, 0.00000000, -0.14142136],
                    [0.70710678, 0.00000000, 0.56568542],
                    [-0.70710678, 0.00000000, 0.56568542],
                ],
                [
                    [0.00000000, 0.00000000, -0.06392934],
                    [0.76924532, 0.00000000, 0.52693942],
                    [-0.76924532, 0.00000000, 0.52693942],
                ],
            ]),
            "angstrom"
        )
        self.assertTrue(np.array_equal(np.array([positions[0], positions[-1]]), expected_pos))

        # Test labels
        scc_indices = self.results["frame_sequence_local_frames_ref"]
        sccs = self.results["section_single_configuration_calculation"]
        systems = self.results["section_system"]
        labels = []
        for index in scc_indices:
            scc = sccs[index]
            system_ref = scc["single_configuration_calculation_to_system_ref"]
            system = systems[system_ref]
            i_labels = system["atom_labels"]
            labels.append(i_labels)
        expected_labels = np.array(4*["O", "H", "H"]).reshape(4, 3)
        self.assertTrue(np.array_equal(labels, expected_labels))

    def test_sampling_method(self):
        result = self.results["sampling_method"]
        self.assertEqual(result, "geometry_optimization")

    def test_geometry_optimization_threshold_force(self):
        result = self.results["geometry_optimization_threshold_force"]
        expected_result = convert_unit(0.000450, "bohr^-1*hartree")
        self.assertEqual(result, expected_result)

    def test_geometry_optimization_energy_change(self):
        result = self.results["geometry_optimization_energy_change"]
        expected_result = convert_unit(5.0E-06, "hartree")
        self.assertEqual(result, expected_result)

    def test_geometry_optimization_geometry_change(self):
        result = self.results["geometry_optimization_geometry_change"]
        expected_result = convert_unit(0.001800, "bohr")
        self.assertEqual(result, expected_result)

    def test_atom_forces(self):
        result = self.results["atom_forces"]
        expected_start = convert_unit(
            -np.array([
                [0.000000, 0.000000, -0.056081],
                [-0.006520, 0.000000, 0.028040],
                [0.006520, 0.000000, 0.028040],
            ]),
            "hartree/bohr"
        )
        self.assertTrue(np.array_equal(result[0, :, :], expected_start))
        expected_end = convert_unit(
            -np.array([
                [0.000000, 0.000000, -0.000037],
                [0.000006, 0.000000, 0.000018],
                [-0.000006, 0.000000, 0.000018],
            ]),
            "hartree/bohr"
        )
        self.assertTrue(np.array_equal(result[-1, :, :], expected_end))
        self.assertEqual(len(result), 4)

    def test_energy_total(self):
        result = self.results["energy_total"]
        np.set_printoptions(precision=12)
        expected_start = convert_unit(-76.429418611381, "hartree")
        self.assertTrue(np.array_equal(result[0], expected_start))

    def test_frame_sequence_to_sampling_ref(self):
        result = self.results["frame_sequence_to_sampling_ref"]
        self.assertEqual(result, 0)

    def test_frame_sequence_local_frames_ref(self):
        result = self.results["frame_sequence_local_frames_ref"]
        expected_result = np.array([0, 2, 4, 6])
        self.assertTrue(np.array_equal(result, expected_result))


#===============================================================================
class TestDFTGaussianMD(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.results = get_results("dft_gaussian/md", "section_run")
        cls.time = convert_unit(
            np.array([
                0.241888,
                0.483777,
                0.725665,
                0.967554,
                1.209442,
            ]),
            "fs"
        )
        cls.kin = convert_unit(
            np.array([
                0.001588,
                0.001316,
                0.001306,
                0.000875,
                0.000360,
            ]),
            "hartree"
        )
        cls.pot = convert_unit(
            np.array([
                -76.325047,
                -76.324974,
                -76.324885,
                -76.324800,
                -76.324755,
            ]),
            "hartree"
        )
        cls.cons = convert_unit(
            np.array([
                -76.323459,
                -76.323657,
                -76.323578,
                -76.323925,
                -76.324394,
            ]),
            "hartree"
        )
        cls.temp = convert_unit(
            np.array([
                334.35,
                277.07,
                275.04,
                184.29,
                75.89,
            ]),
            "K"
        )

    def get_system(self, index):
        scc = self.get_scc(index)
        system_ref = scc["single_configuration_calculation_to_system_ref"]
        system = self.results["section_system"][system_ref]
        return system

    def get_scc(self, index):
        sample_refs = self.results["frame_sequence_local_frames_ref"]
        sccs = self.results["section_single_configuration_calculation"]
        scc = sccs[sample_refs[index]]
        return scc

    def test_configuration_periodic_dimensions(self):
        result = self.results["configuration_periodic_dimensions"][0]
        self.assertTrue(np.array_equal(result, np.array([False, False, False])))

    def test_single_configuration_to_calculation_method(self):
        result = self.results["single_configuration_to_calculation_method_ref"]
        self.assertTrue(np.array_equal(result, np.array(6*[0])))

    def test_electronic_structure_method(self):
        result = self.results["electronic_structure_method"]
        self.assertEqual(result, "DFT")

    def test_xc_functional(self):
        result = self.results["XC_functional"]
        self.assertEqual(result, "1.0*HYB_GGA_XC_PBEH")

    def test_sampling_method(self):
        result = self.results["sampling_method"]
        self.assertEqual(result, "molecular_dynamics")

    def test_ensemble_type(self):
        result = self.results["ensemble_type"]
        self.assertEqual(result, "NVT")

    def test_frame_sequence_to_sampling_ref(self):
        result = self.results["frame_sequence_to_sampling_ref"]
        self.assertEqual(result, 0)

    def test_frame_sequence_local_frames_ref(self):
        result = self.results["frame_sequence_local_frames_ref"]
        self.assertTrue(np.array_equal(result, np.array(range(1, 6))))

    def test_number_of_frames_in_sequence(self):
        result = self.results["number_of_frames_in_sequence"]
        self.assertEqual(result, 5)

    def test_frame_sequence_potential_energy(self):
        result = self.results["frame_sequence_potential_energy"]
        self.assertTrue(np.array_equal(result, self.pot))

    def test_frame_sequence_kinetic_energy(self):
        result = self.results["frame_sequence_kinetic_energy"]
        self.assertTrue(np.array_equal(result, self.kin))

    def test_frame_sequence_temperature(self):
        result = self.results["frame_sequence_temperature"]
        self.assertTrue(np.array_equal(result, self.temp))

    def test_frame_sequence_time(self):
        result = self.results["frame_sequence_time"]
        self.assertTrue(np.array_equal(result, self.time))

    def test_frame_sequence_potential_energy_stats(self):
        result = self.results["frame_sequence_potential_energy_stats"]
        expected_result = np.array([self.pot.mean(), self.pot.std()])
        self.assertTrue(np.allclose(result[0], expected_result[0], rtol=0, atol=0.00001e-18))
        self.assertTrue(np.allclose(result[1], expected_result[1], rtol=0, atol=0.00001e-20))

    def test_frame_sequence_kinetic_energy_stats(self):
        result = self.results["frame_sequence_kinetic_energy_stats"]
        expected_result = np.array([self.kin.mean(), self.kin.std()])
        self.assertTrue(np.allclose(result[0], expected_result[0], rtol=0, atol=0.00001e-18))
        self.assertTrue(np.allclose(result[1], expected_result[1], rtol=0, atol=0.00001e-20))

    def test_frame_sequence_temperature_stats(self):
        result = self.results["frame_sequence_temperature_stats"]
        expected_result = np.array([self.temp.mean(), self.temp.std()])
        self.assertTrue(np.allclose(result[0], expected_result[0], rtol=0, atol=0.001))
        self.assertTrue(np.allclose(result[1], expected_result[1], rtol=0, atol=0.0001))

    def test_atom_positions(self):
        first_system = self.get_system(0)
        first_pos = first_system["atom_positions"]
        last_system = self.get_system(-1)
        last_pos = last_system["atom_positions"]
        expected_start = convert_unit(
            np.array([
                [-0.000000,  -0.032407,   0.213730],
                [ 0.000000,   1.547303,  -0.646826],
                [-0.000000,  -1.283238,  -1.058258],
            ]),
            "bohr"
        )
        expected_end = convert_unit(
            np.array([
                [-0.000000,  -0.034144,   0.212417],
                [ 0.000000,   1.583695,  -0.644729],
                [-0.000000,  -1.292061,  -1.039511],
            ]),
            "bohr"
        )

        self.assertTrue(np.array_equal(first_pos, expected_start))
        self.assertTrue(np.array_equal(last_pos, expected_end))

    def test_atom_forces(self):
        first_force = self.get_scc(0)["atom_forces"]
        last_force = self.get_scc(-1)["atom_forces"]
        expected_start = convert_unit(
            -np.array([
                [ 0.000000,  -0.003686,  -0.024792],
                [-0.000000,  -0.009261,   0.007954],
                [ 0.000000,   0.012947,   0.016838],
            ]),
            "forceAu"
        )
        expected_end = convert_unit(
            -np.array([
                [-0.000000,  -0.023297,  -0.023732],
                [-0.000000,   0.008095,   0.001352],
                [ 0.000000,   0.015202,   0.022380],
            ]),
            "forceAu"
        )

        self.assertTrue(np.array_equal(first_force, expected_start))
        self.assertTrue(np.array_equal(last_force, expected_end))

    # def test_atom_velocities(self):
        # result = self.results["atom_velocities"]
        # expected_start = convert_unit(
            # np.array([
                # [0.00039772295627, 0.00024115257177, 0.00026132422738],
                # [-0.00039772295627, -0.00024115257177, -0.00026132422738],
            # ]),
            # "bohr/(hbar/hartree)"
        # )
        # expected_end = convert_unit(
            # np.array([
                # [0.00094644268934, 0.00023563385430, 0.00025534388718],
                # [-0.00094644268934, -0.00023563385430, -0.00025534388718],
            # ]),
            # "bohr/(hbar/hartree)"
        # )

        # self.assertTrue(np.array_equal(result[0, :], expected_start))
        # self.assertTrue(np.array_equal(result[-1, :], expected_end))


#===============================================================================
class TestDFTGaussianXCFunctional(unittest.TestCase):
    """Tests that the XC functionals can be properly parsed.
    """
    def test_blyp(self):
        xc = get_result("dft_gaussian/functionals/blyp", "XC_functional")
        self.assertEqual(xc, "1.0*GGA_C_LYP+1.0*GGA_X_B88")

    def test_b3lyp(self):
        xc = get_result("dft_gaussian/functionals/b3lyp", "XC_functional")
        self.assertEqual(xc, "1.0*HYB_GGA_XC_B3LYP")

    def test_pbe(self):
        xc = get_result("dft_gaussian/functionals/pbe", "XC_functional")
        self.assertEqual(xc, "1.0*GGA_C_PBE+1.0*GGA_X_PBE")

    def test_pbe0(self):
        xc = get_result("dft_gaussian/functionals/pbe0", "XC_functional")
        self.assertEqual(xc, "1.0*HYB_GGA_XC_PBEH")

    def test_bp86(self):
        xc = get_result("dft_gaussian/functionals/bp86", "XC_functional")
        self.assertEqual(xc, "1.0*GGA_C_P86+1.0*GGA_X_B88")

    def test_bp91(self):
        xc = get_result("dft_gaussian/functionals/bp91", "XC_functional")
        self.assertEqual(xc, "1.0*GGA_C_PW91+1.0*GGA_X_B88")

    def test_pw91(self):
        xc = get_result("dft_gaussian/functionals/pw91", "XC_functional")
        self.assertEqual(xc, "1.0*GGA_C_PW91+1.0*GGA_X_PW91")

    def test_bechehandh(self):
        xc = get_result("dft_gaussian/functionals/beckehandh", "XC_functional")
        self.assertEqual(xc, "1.0*HYB_GGA_XC_BHANDH")

    def test_olyp(self):
        xc = get_result("dft_gaussian/functionals/olyp", "XC_functional")
        self.assertEqual(xc, "1.0*GGA_C_LYP+1.0*GGA_X_OPTX")

    def test_hcth120(self):
        xc = get_result("dft_gaussian/functionals/hcth120", "XC_functional")
        self.assertEqual(xc, "1.0*GGA_XC_HCTH_120")

    def test_hcth147(self):
        xc = get_result("dft_gaussian/functionals/hcth147", "XC_functional")
        self.assertEqual(xc, "1.0*GGA_XC_HCTH_147")

    def test_hcth407(self):
        xc = get_result("dft_gaussian/functionals/hcth407", "XC_functional")
        self.assertEqual(xc, "1.0*GGA_XC_HCTH_407")

    def test_tpss(self):
        xc = get_result("dft_gaussian/functionals/tpss", "XC_functional")
        self.assertEqual(xc, "1.0*MGGA_C_TPSS+1.0*MGGA_X_TPSS")


#===============================================================================
class TestDFTPWEnergy(unittest.TestCase):
    """Tests that the parser can handle plane-wave DFT energy calculations.
    """

    @classmethod
    def setUpClass(cls):
        cls.results = get_results("dft_pw/energy", "section_run")
        # cls.results.print_summary()

    def test_program_name(self):
        result = self.results["program_name"]
        self.assertEqual(result, "NWChem")

    def test_configuration_periodic_dimensions(self):
        result = self.results["configuration_periodic_dimensions"][0]
        self.assertTrue(np.array_equal(result, np.array([True, True, True])))

    def test_program_version(self):
        result = self.results["program_version"]
        self.assertEqual(result, "6.6")

    def test_total_charge(self):
        charge = self.results["total_charge"][0]
        self.assertEqual(charge, 0)

    def test_electronic_structure_method(self):
        result = self.results["electronic_structure_method"][0]
        self.assertEqual(result, "DFT")

    def test_simulation_cell(self):
        result = self.results["simulation_cell"][0]
        self.assertTrue(np.array_equal(result, convert_unit(np.array(
            [
                [20.0, 0.0, 0.0],
                [0.0, 20.0, 0.0],
                [0.0, 0.0, 20.0],
            ]), "angstrom")
        ))

    # def test_atom_labels(self):
        # atom_labels = self.results["atom_labels"]
        # expected_labels = np.array(["O", "H", "H"])
        # self.assertTrue(np.array_equal(atom_labels, expected_labels))

    # def test_atom_positions(self):
        # atom_position = self.results["atom_positions"]
        # expected_position = convert_unit(np.array(
            # [
                # [0.00000000, 0.00000000, -0.11817375],
                # [0.76924532, 0.00000000, 0.47269501],
                # [-0.76924532, 0.00000000, 0.47269501],
            # ]
        # ), "angstrom")
        # self.assertTrue(np.array_equal(atom_position, expected_position))

    # def test_number_of_atoms(self):
        # n_atoms = self.results["number_of_atoms"]
        # self.assertEqual(n_atoms, 3)

    # def test_energy_total(self):
        # result = self.results["energy_total"]
        # expected_result = convert_unit(np.array(-76.436222730188), "hartree")
        # self.assertTrue(np.array_equal(result, expected_result))

    # def test_energy_x(self):
        # result = self.results["energy_X"]
        # expected_result = convert_unit(np.array(-9.025345841743), "hartree")
        # self.assertTrue(np.array_equal(result, expected_result))

    # def test_energy_c(self):
        # result = self.results["energy_C"]
        # expected_result = convert_unit(np.array(-0.328011552453), "hartree")
        # self.assertTrue(np.array_equal(result, expected_result))

    # def test_energy_total_scf_iteration(self):
        # result = self.results["energy_total_scf_iteration"]
        # # Test the first and last energies
        # expected_result = convert_unit(np.array(
            # [
                # [-76.3916403957],
                # [-76.4362227302],
            # ]), "hartree")
        # self.assertTrue(np.array_equal(np.array([[result[0]], [result[-1]]]), expected_result))

    # def test_energy_change_scf_iteration(self):
        # result = self.results["energy_change_scf_iteration"]
        # expected_result = convert_unit(np.array(
            # [
                # [-8.55E+01],
                # [-3.82E-07],
            # ]), "hartree")
        # self.assertTrue(np.array_equal(np.array([[result[0]], [result[-1]]]), expected_result))

    # def test_scf_max_iteration(self):
        # result = self.results["scf_max_iteration"]
        # self.assertEqual(result, 50)

    # def test_scf_threshold_energy_change(self):
        # result = self.results["scf_threshold_energy_change"]
        # self.assertEqual(result, convert_unit(1.00E-06, "hartree"))

    # def test_scf_dft_number_of_iterations(self):
        # result = self.results["number_of_scf_iterations"]
        # self.assertEqual(result, 6)

    # def test_spin_target_multiplicity(self):
        # multiplicity = self.results["spin_target_multiplicity"]
        # self.assertEqual(multiplicity, 1)

    # def test_single_configuration_to_calculation_method_ref(self):
        # result = self.results["single_configuration_to_calculation_method_ref"]
        # self.assertEqual(result, 0)

    # def test_single_configuration_calculation_to_system_description_ref(self):
        # result = self.results["single_configuration_calculation_to_system_ref"]
        # self.assertEqual(result, 0)


#===============================================================================
if __name__ == '__main__':
    suites = []
    suites.append(unittest.TestLoader().loadTestsFromTestCase(TestDFTGaussianEnergy))
    suites.append(unittest.TestLoader().loadTestsFromTestCase(TestDFTGaussianForce))
    suites.append(unittest.TestLoader().loadTestsFromTestCase(TestDFTGaussianGeoOpt))
    suites.append(unittest.TestLoader().loadTestsFromTestCase(TestDFTGaussianXCFunctional))
    suites.append(unittest.TestLoader().loadTestsFromTestCase(TestDFTGaussianMD))

    suites.append(unittest.TestLoader().loadTestsFromTestCase(TestDFTPWEnergy))

    alltests = unittest.TestSuite(suites)
    unittest.TextTestRunner(verbosity=0).run(alltests)
