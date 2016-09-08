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
class TestDFTEnergy(unittest.TestCase):
    """Tests that the parser can handle DFT energy calculations.
    """

    @classmethod
    def setUpClass(cls):
        cls.results = get_results("dft/energy", "section_run")
        # cls.results.print_summary()

    def test_program_name(self):
        result = self.results["program_name"]
        self.assertEqual(result, "NWChem")

    def test_program_version(self):
        result = self.results["program_version"]
        self.assertEqual(result, "6.6")

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

    # def test_stress_tensor(self):
        # result = self.results["stress_tensor"]
        # expected_result = convert_unit(
            # np.array([
                # [7.77640934, -0.00000098, -0.00000099],
                # [-0.00000098, 7.77640935, -0.00000101],
                # [-0.00000099, -0.00000101, 7.77640935],
            # ]),
            # "GPa"
        # )
        # self.assertTrue(np.array_equal(result, expected_result))

    # def test_stress_tensor_eigenvalues(self):
        # result = self.results["x_cp2k_stress_tensor_eigenvalues"]
        # expected_result = convert_unit(np.array([7.77640735, 7.77641033, 7.77641036]), "GPa")
        # self.assertTrue(np.array_equal(result, expected_result))

    # def test_stress_tensor_eigenvectors(self):
        # result = self.results["x_cp2k_stress_tensor_eigenvectors"]
        # expected_result = np.array([
            # [0.57490332, -0.79965737, -0.17330395],
            # [0.57753686, 0.54662171, -0.60634634],
            # [0.57960102, 0.24850110, 0.77608624],
        # ])
        # self.assertTrue(np.array_equal(result, expected_result))

    # def test_stress_tensor_determinant(self):
        # result = self.results["x_cp2k_stress_tensor_determinant"]
        # expected_result = convert_unit(4.70259243E+02, "GPa^3")
        # self.assertTrue(np.array_equal(result, expected_result))

    # def test_stress_tensor_one_third_of_trace(self):
        # result = self.results["x_cp2k_stress_tensor_one_third_of_trace"]
        # expected_result = convert_unit(7.77640934E+00, "GPa")
        # self.assertTrue(np.array_equal(result, expected_result))


#===============================================================================
class TestDFTForce(unittest.TestCase):
    """Tests that the parser can handle DFT force calculations.
    """
    @classmethod
    def setUpClass(cls):
        cls.results = get_results("dft/force", "section_run")

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
class TestDFTGeoOpt(unittest.TestCase):
    """Tests that the parser can handle DFT geometry optimizations.
    """
    @classmethod
    def setUpClass(cls):
        cls.results = get_results("dft/geo_opt", "section_run")

    def test_frame_sequence(self):
        sequence = self.results["section_frame_sequence"][0]

        # Number of frames
        n_frames = self.results["number_of_frames_in_sequence"]
        self.assertEqual(n_frames, 4)

        # Potential energy
        pot_ener = sequence["frame_sequence_potential_energy"]
        # print(pot_ener)
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
        labels = self.results["atom_labels"]
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
        self.assertEqual(len(result), 4)
        self.assertTrue(np.array_equal(result[-1, :, :], expected_end))

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
        expected_result = np.array([0, 1, 2, 3])
        self.assertTrue(np.array_equal(result, expected_result))


#===============================================================================
# class TestMD(unittest.TestCase):
    # @classmethod
    # def setUpClass(cls):
        # cls.results = get_results("md/nve", "section_run")
        # cls.temp = convert_unit(
            # np.array([
                # 110.096,
                # 232.496,
                # 351.956,
                # 412.578,
                # 393.180,
            # ]),
            # "K"
        # )
        # cls.cons = convert_unit(
            # np.array([
                # -1.0970967730,
                # -1.0975238350,
                # -1.0977293448,
                # -1.0977368045,
                # -1.0975921059,
            # ]),
            # "hartree"
        # )
        # cls.pot = convert_unit(
            # np.array([
                # -1.1023492072,
                # -1.1128688938,
                # -1.1216882365,
                # -1.1256188624,
                # -1.1245335482,
            # ]),
            # "hartree"
        # )
        # cls.kin = convert_unit(
            # np.array([
                # -1.1018262261,
                # -1.1117644858,
                # -1.1200163669,
                # -1.1236590243,
                # -1.1226658551,
            # ]) -
            # np.array([
                # -1.1023492072,
                # -1.1128688938,
                # -1.1216882365,
                # -1.1256188624,
                # -1.1245335482,
            # ]),
            # "hartree"
        # )

    # def test_number_of_atoms(self):
        # result = self.results["number_of_atoms"]
        # expected_result = np.array(5*[2])
        # self.assertTrue(np.array_equal(result, expected_result))

    # def test_ensemble_type(self):
        # result = self.results["ensemble_type"]
        # self.assertEqual(result, "NVE")

    # def test_sampling_method(self):
        # result = self.results["sampling_method"]
        # self.assertEqual(result, "molecular_dynamics")

    # def test_number_of_frames_in_sequence(self):
        # result = self.results["number_of_frames_in_sequence"]
        # self.assertEqual(result, 5)

    # def test_atom_positions(self):
        # result = self.results["atom_positions"]
        # expected_start = convert_unit(
            # np.array([
                # [0.371489, 0.000511, 0.000554],
                # [-0.371489, -0.000511, -0.000554],
            # ]),
            # "angstrom"
        # )
        # expected_end = convert_unit(
            # np.array([
                # [0.378523, 0.002532, 0.002744],
                # [-0.378523, -0.002532, -0.002744],
            # ]),
            # "angstrom"
        # )
        # self.assertTrue(np.array_equal(result[0, :], expected_start))
        # self.assertTrue(np.array_equal(result[-1, :], expected_end))

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

    # def test_atom_forces(self):
        # result = self.results["atom_forces"]
        # expected_start = convert_unit(
            # np.array([
                # [0.15293653991241, -0.00036559789218, -0.00039617900820],
                # [-0.15293653991238, 0.00036559789217, 0.00039617900820],
            # ]),
            # "forceAu"
        # )
        # expected_end = convert_unit(
            # np.array([
                # [-0.03595092814462, -0.00079338139843, -0.00085974499854],
                # [0.03595092814462, 0.00079338139843, 0.00085974499854],
            # ]),
            # "forceAu"
        # )

        # self.assertTrue(np.array_equal(result[0, :], expected_start))
        # self.assertTrue(np.array_equal(result[-1, :], expected_end))

    # def test_frame_sequence_potential_energy(self):
        # result = self.results["frame_sequence_potential_energy"]
        # self.assertTrue(np.array_equal(result, self.pot))

    # def test_frame_sequence_kinetic_energy(self):
        # result = self.results["frame_sequence_kinetic_energy"]
        # self.assertTrue(np.array_equal(result, self.kin))

    # def test_frame_sequence_conserved_quantity(self):
        # result = self.results["frame_sequence_conserved_quantity"]
        # self.assertTrue(np.array_equal(result, self.cons))

    # def test_frame_sequence_temperature(self):
        # result = self.results["frame_sequence_temperature"]
        # self.assertTrue(np.array_equal(result, self.temp))

    # def test_frame_sequence_time(self):
        # result = self.results["frame_sequence_time"]
        # expected_result = convert_unit(
            # np.array([
                # 4,
                # 8,
                # 12,
                # 16,
                # 20,
            # ]),
            # "hbar/hartree"
        # )
        # self.assertTrue(np.array_equal(result, expected_result))

    # def test_frame_sequence_potential_energy_stats(self):
        # result = self.results["frame_sequence_potential_energy_stats"]
        # expected_result = np.array([self.pot.mean(), self.pot.std()])
        # self.assertTrue(np.allclose(result[0], expected_result[0], rtol=0, atol=0.00001e-18))
        # self.assertTrue(np.allclose(result[1], expected_result[1], rtol=0, atol=0.00001e-20))

    # def test_frame_sequence_kinetic_energy_stats(self):
        # result = self.results["frame_sequence_kinetic_energy_stats"]
        # expected_result = np.array([self.kin.mean(), self.kin.std()])
        # self.assertTrue(np.allclose(result[0], expected_result[0], rtol=0, atol=0.00001e-18))
        # self.assertTrue(np.allclose(result[1], expected_result[1], rtol=0, atol=0.00001e-20))

    # def test_frame_sequence_conserved_quantity_stats(self):
        # result = self.results["frame_sequence_conserved_quantity_stats"]
        # expected_result = np.array([self.cons.mean(), self.cons.std()])
        # self.assertTrue(np.allclose(result[0], expected_result[0], rtol=0, atol=0.00001e-18))
        # self.assertTrue(np.allclose(result[1], expected_result[1], rtol=0, atol=0.00001e-20))

    # def test_frame_sequence_temperature_stats(self):
        # result = self.results["frame_sequence_temperature_stats"]
        # expected_result = np.array([self.temp.mean(), self.temp.std()])
        # self.assertTrue(np.allclose(result[0], expected_result[0], rtol=0, atol=0.001))
        # self.assertTrue(np.allclose(result[1], expected_result[1], rtol=0, atol=0.0001))


#===============================================================================
# class TestMDTrajFormats(unittest.TestCase):

    # def test_dcd(self):
        # results = get_results("md/dcd", "section_run")
        # positions = results["atom_positions"]
        # expected_start = convert_unit(
            # np.array([
                # [0.70201340784773, 0.00096620207776, 0.00104702184853],
                # [-0.70201340784773, -0.00096620207776, -0.00104702184853],
            # ]),
            # "bohr"
        # )
        # expected_end = convert_unit(
            # np.array([
                # [0.71530475161473, 0.00478536889215, 0.00518564998016],
                # [-0.71530475161473, -0.00478536889215, -0.00518564998016],
            # ]),
            # "bohr"
        # )
        # self.assertTrue(np.allclose(positions[0, :], expected_start, rtol=0, atol=0.0000001e-11))
        # self.assertTrue(np.allclose(positions[-1, :], expected_end, rtol=0, atol=0.0000001e-11))

    # def test_xyz(self):
        # results = get_results("md/xyz", "section_run")
        # positions = results["atom_positions"]
        # expected_start = convert_unit(
            # np.array([
                # [0.70201340784773, 0.00096620207776, 0.00104702184853],
                # [-0.70201340784773, -0.00096620207776, -0.00104702184853],
            # ]),
            # "bohr"
        # )
        # expected_end = convert_unit(
            # np.array([
                # [0.71530475161473, 0.00478536889215, 0.00518564998016],
                # [-0.71530475161473, -0.00478536889215, -0.00518564998016],
            # ]),
            # "bohr"
        # )
        # self.assertTrue(np.allclose(positions[0, :], expected_start, rtol=0, atol=0.00001e-11))
        # self.assertTrue(np.allclose(positions[-1, :], expected_end, rtol=0, atol=0.00001e-11))

    # def test_trajectory(self):
        # results = get_results("md/trajectory", "section_run")
        # positions = results["atom_positions"]
        # expected_start = convert_unit(
            # np.array([
                # [0.70201340784773, 0.00096620207776, 0.00104702184853],
                # [-0.70201340784773, -0.00096620207776, -0.00104702184853],
            # ]),
            # "bohr"
        # )
        # expected_end = convert_unit(
            # np.array([
                # [0.71530475161473, 0.00478536889215, 0.00518564998016],
                # [-0.71530475161473, -0.00478536889215, -0.00518564998016],
            # ]),
            # "bohr"
        # )
        # self.assertTrue(np.allclose(positions[0, :], expected_start, rtol=0, atol=0.00001e-11))
        # self.assertTrue(np.allclose(positions[-1, :], expected_end, rtol=0, atol=0.00001e-11))

    # def test_ftrajectory(self):
        # results = get_results("md/ftrajectory", "section_run")
        # positions = results["atom_positions"]
        # expected_start = convert_unit(
            # np.array([
                # [0.70201340784773, 0.00096620207776, 0.00104702184853],
                # [-0.70201340784773, -0.00096620207776, -0.00104702184853],
            # ]),
            # "bohr"
        # )
        # expected_end = convert_unit(
            # np.array([
                # [0.71530475161473, 0.00478536889215, 0.00518564998016],
                # [-0.71530475161473, -0.00478536889215, -0.00518564998016],
            # ]),
            # "bohr"
        # )
        # self.assertTrue(np.allclose(positions[0, :], expected_start, rtol=0, atol=0.00001e-11))
        # self.assertTrue(np.allclose(positions[-1, :], expected_end, rtol=0, atol=0.00001e-11))


#===============================================================================
class TestXCFunctional(unittest.TestCase):
    """Tests that the XC functionals can be properly parsed.
    """

    # def test_lda(self):
        # xc = get_result("xc_functional/lda", "XC_functional")
        # self.assertEqual(xc, "1*LDA_XC_TETER93")

    # def test_blyp(self):
        # xc = get_result("xc_functional/blyp", "XC_functional")
        # self.assertEqual(xc, "1*GGA_C_LYP+1*GGA_X_B88")

    def test_b3lyp(self):
        xc = get_result("dft/functionals/b3lyp", "XC_functional")
        self.assertEqual(xc, "1.0*HYB_GGA_XC_B3LYP")

    def test_pbe(self):
        xc = get_result("dft/functionals/pbe", "XC_functional")
        self.assertEqual(xc, "1.0*GGA_C_PBE+1.0*GGA_X_PBE")

    # def test_olyp(self):
        # xc = get_result("xc_functional/olyp", "XC_functional")
        # self.assertEqual(xc, "1*GGA_C_LYP+1*GGA_X_OPTX")

    # def test_hcth(self):
        # xc = get_result("xc_functional/hcth", "XC_functional")
        # self.assertEqual(xc, "1*GGA_XC_HCTH_120")

    def test_pbe0(self):
        xc = get_result("dft/functionals/pbe0", "XC_functional")
        self.assertEqual(xc, "1.0*HYB_GGA_XC_PBEH")

    # def test_bp(self):
        # xc = get_result("xc_functional/bp", "XC_functional")
        # self.assertEqual(xc, "1*GGA_C_P86+1*GGA_X_B88")

    # def test_xlyp(self):
        # xc = get_result("xc_functional/xlyp", "XC_functional")
        # self.assertEqual(xc, "1*GGA_XC_XLYP")

    # def test_pbes(self):
        # xc = get_result("xc_functional/pbes", "XC_functional")
        # self.assertEqual(xc, "1*GGA_C_PBE_SOL+1*GGA_X_PBE_SOL")

    # def test_revpbe(self):
        # xc = get_result("xc_functional/revpbe", "XC_functional")
        # self.assertEqual(xc, "1*GGA_C_PBE+1*GGA_X_PBE_R")

    # def test_tpss(self):
        # xc = get_result("xc_functional/tpss", "XC_functional")
        # self.assertEqual(xc, "1*MGGA_C_TPSS+1*MGGA_X_TPSS")

    # def test_b1lyp(self):
        # xc = get_result("xc_functional/b1lyp", "XC_functional")
        # self.assertEqual(xc, "1*HYB_GGA_XC_B1LYP")

    # def test_x3lyp(self):
        # xc = get_result("xc_functional/x3lyp", "XC_functional")
        # self.assertEqual(xc, "1*HYB_GGA_XC_X3LYP")

    # def test_hse06(self):
        # xc = get_result("xc_functional/hse06", "XC_functional")
        # self.assertEqual(xc, "1*HYB_GGA_XC_HSE06")


# #===============================================================================
# class TestErrors(unittest.TestCase):
    # """Test misc. error stuations which may occur during the parsing.
    # """
    # def test_no_file(self):
        # self.assertRaises(IOError, get_result, "errors/no_file", "XC_functional")

    # def test_invalid_file(self):
        # self.assertRaises(RuntimeError, get_result, "errors/invalid_file", "XC_functional")

    # def test_invalid_run_type(self):
        # self.assertRaises(KeyError, get_result, "errors/invalid_run_type", "XC_functional")

    # def test_unknown_version(self):
        # get_result("errors/unknown_version", "XC_functional")

    # def test_unknown_input_keyword(self):
        # get_result("errors/unknown_input_keyword", "XC_functional")

    # def test_unknown_input_section(self):
        # get_result("errors/unknown_input_section", "XC_functional")

    # def test_unknown_input_section_parameter(self):
        # get_result("errors/unknown_input_section_parameter", "XC_functional")


# #===============================================================================
# class TestSCFConvergence(unittest.TestCase):
    # """Tests whether the convergence status and number of SCF step can be
    # parsed correctly.
    # """

    # def test_converged(self):
        # result = get_result("convergence/converged", "single_configuration_calculation_converged")
        # self.assertTrue(result)

    # def test_non_converged(self):
        # result = get_result("convergence/non_converged", "single_configuration_calculation_converged")
        # self.assertFalse(result)


#===============================================================================
# class TestSelfInteractionCorrectionMethod(unittest.TestCase):
    # """Tests that the self-interaction correction can be properly parsed.
    # """

    # def test_no(self):
        # sic = get_result("sic/no", "self_interaction_correction_method")
        # self.assertEqual(sic, "")

    # def test_ad(self):
        # sic = get_result("sic/ad", "self_interaction_correction_method")
        # self.assertEqual(sic, "SIC_AD")

    # def test_explicit_orbitals(self):
        # sic = get_result("sic/explicit_orbitals", "self_interaction_correction_method")
        # self.assertEqual(sic, "SIC_EXPLICIT_ORBITALS")

    # def test_mauri_spz(self):
        # sic = get_result("sic/mauri_spz", "self_interaction_correction_method")
        # self.assertEqual(sic, "SIC_MAURI_SPZ")

    # def test_mauri_us(self):
        # sic = get_result("sic/mauri_us", "self_interaction_correction_method")
        # self.assertEqual(sic, "SIC_MAURI_US")


# #===============================================================================
# class TestStressTensorMethods(unittest.TestCase):
    # """Tests that the stress tensor can be properly parsed for different
    # calculation methods.
    # """
    # def test_none(self):
        # get_results("stress_tensor/none", "section_stress_tensor")

    # def test_analytical(self):
        # results = get_results("stress_tensor/analytical", ["stress_tensor_method", "stress_tensor"])
        # method = results["stress_tensor_method"]
        # results["stress_tensor"]
        # self.assertEqual(method, "Analytical")

    # def test_numerical(self):
        # results = get_results("stress_tensor/numerical", ["stress_tensor_method", "stress_tensor"])
        # method = results["stress_tensor_method"]
        # results["stress_tensor"]
        # self.assertEqual(method, "Numerical")

    # def test_diagonal_analytical(self):
        # results = get_results("stress_tensor/diagonal_analytical", ["stress_tensor_method", "stress_tensor"])
        # method = results["stress_tensor_method"]
        # results["stress_tensor"]
        # self.assertEqual(method, "Diagonal analytical")

    # def test_diagonal_numerical(self):
        # results = get_results("stress_tensor/diagonal_numerical", ["stress_tensor_method", "stress_tensor"])
        # method = results["stress_tensor_method"]
        # results["stress_tensor"]
        # self.assertEqual(method, "Diagonal numerical")


# # ===============================================================================
# class TestGeoOptTrajFormats(unittest.TestCase):

    # def test_xyz(self):

        # result = get_result("geo_opt/geometry_formats/xyz", "atom_positions", optimize=True)
        # expected_start = convert_unit(
            # np.array([
                # [12.2353220000, 1.3766420000, 10.8698800000],
                # [12.4175624065, 2.2362390825, 11.2616392180],
                # [11.9271777126, 1.5723402996, 10.0115089094],
            # ]),
            # "angstrom"
        # )
        # expected_end = convert_unit(
            # np.array([
                # [12.2353220000, 1.3766420000, 10.8698800000],
                # [12.4957995882, 2.2307218433, 11.3354453867],
                # [11.9975764125, 1.5747996320, 10.0062529540],
            # ]),
            # "angstrom"
        # )
        # result_start = result[0,:,:]
        # result_end = result[-1,:,:]
        # self.assertTrue(np.array_equal(result_start, expected_start))
        # self.assertTrue(np.array_equal(result_end, expected_end))

    # def test_pdb(self):
        # result = get_result("geo_opt/geometry_formats/pdb", "atom_positions", optimize=True)
        # expected_start = convert_unit(
            # np.array([
                # [12.235, 1.377, 10.870],
                # [12.418, 2.236, 11.262],
                # [11.927, 1.572, 10.012],
            # ]),
            # "angstrom"
        # )
        # expected_end = convert_unit(
            # np.array([
                # [12.235, 1.377, 10.870],
                # [12.496, 2.231, 11.335],
                # [11.998, 1.575, 10.006],
            # ]),
            # "angstrom"
        # )
        # result_start = result[0,:,:]
        # result_end = result[-1,:,:]
        # self.assertTrue(np.array_equal(result_start, expected_start))
        # self.assertTrue(np.array_equal(result_end, expected_end))

    # def test_dcd(self):
        # result = get_result("geo_opt/geometry_formats/dcd", "atom_positions", optimize=True)
        # frames = result.shape[0]
        # self.assertEqual(frames, 7)


# #===============================================================================
# class TestGeoOptOptimizers(unittest.TestCase):

    # def test_bfgs(self):
        # result = get_result("geo_opt/bfgs", "geometry_optimization_method")
        # self.assertEqual(result, "bfgs")

    # def test_lbfgs(self):
        # result = get_result("geo_opt/lbfgs", "geometry_optimization_method")
        # self.assertEqual(result, "bfgs")


# #===============================================================================
# class TestGeoOptTrajectory(unittest.TestCase):

    # def test_each_and_add_last(self):
        # """Test that the EACH and ADD_LAST settings affect the parsing
        # correctly.
        # """
        # results = get_results("geo_opt/each")

        # single_conf = results["section_single_configuration_calculation"]
        # systems = results["section_system"]

        # i_conf = 0
        # for calc in single_conf.values():
            # system_index = calc["single_configuration_calculation_to_system_ref"][0]
            # system = systems[system_index]
            # pos = system["atom_positions"]

            # if i_conf == 0 or i_conf == 2 or i_conf == 4:
                # self.assertEqual(pos, None)
            # else:
                # pos = system["atom_positions"][0]
                # if i_conf == 1:
                    # expected_pos = convert_unit(
                        # np.array([
                            # [12.2353220000, 1.3766420000, 10.8698800000],
                            # [12.4618486015, 2.2314871691, 11.3335607388],
                            # [11.9990227122, 1.5776813026, 10.0384213366],
                        # ]),
                        # "angstrom"
                    # )
                    # self.assertTrue(np.array_equal(pos, expected_pos))
                # if i_conf == 3:
                    # expected_pos = convert_unit(
                        # np.array([
                            # [12.2353220000, 1.3766420000, 10.8698800000],
                            # [12.4962705528, 2.2308411983, 11.3355758433],
                            # [11.9975151486, 1.5746309898, 10.0054430868],
                        # ]),
                        # "angstrom"
                    # )
                    # self.assertTrue(np.array_equal(pos, expected_pos))
                # if i_conf == 5:
                    # expected_pos = convert_unit(
                        # np.array([
                            # [12.2353220000, 1.3766420000, 10.8698800000],
                            # [12.4958168364, 2.2307249171, 11.3354322532],
                            # [11.9975556812, 1.5748088251, 10.0062793864],
                        # ]),
                        # "angstrom"
                    # )
                    # self.assertTrue(np.array_equal(pos, expected_pos))

                # if i_conf == 6:
                    # expected_pos = convert_unit(
                        # np.array([
                            # [12.2353220000, 1.3766420000, 10.8698800000],
                            # [12.4958164689, 2.2307248873, 11.3354322515],
                            # [11.9975558616, 1.5748085240, 10.0062792262],
                        # ]),
                        # "angstrom"
                    # )
                    # self.assertTrue(np.array_equal(pos, expected_pos))

            # i_conf += 1


# #===============================================================================
# class TestMDEnsembles(unittest.TestCase):

    # @classmethod
    # def setUpClass(cls):
        # cls.pressure = convert_unit(
            # np.array([
                # -0.192828092559E+04,
                # -0.145371071470E+04,
                # -0.210098903760E+03,
                # 0.167260570313E+04,
                # 0.395562042841E+04,
                # 0.630374855942E+04,
                # 0.836906136786E+04,
                # 0.983216022830E+04,
                # 0.104711540465E+05,
                # 0.102444821550E+05,
                # 0.931695792434E+04,
            # ]),
            # "bar"
        # )

    # def test_nvt(self):
        # results = get_results("md/nvt", "section_run")
        # ensemble = results["ensemble_type"]
        # self.assertEqual(ensemble, "NVT")

    # def test_npt(self):
        # results = get_results("md/npt", "section_run")
        # ensemble = results["ensemble_type"]
        # self.assertEqual(ensemble, "NPT")

        # pressure = results["frame_sequence_pressure"]
        # self.assertTrue(np.array_equal(pressure, self.pressure))

        # pressure_stats = results["frame_sequence_pressure_stats"]
        # expected_pressure_stats = np.array([self.pressure.mean(), self.pressure.std()])
        # self.assertTrue(np.array_equal(pressure_stats, expected_pressure_stats))

        # simulation_cell = results["simulation_cell"]
        # expected_cell_start = convert_unit(
            # np.array(
                # [[
                    # 6.0000000000,
                    # 0.0000000000,
                    # 0.0000000000,
                # ], [
                    # 0.0000000000,
                    # 6.0000000000,
                    # 0.0000000000,
                # ], [
                    # 0.0000000000,
                    # 0.0000000000,
                    # 6.0000000000,
                # ]]),
            # "angstrom"
        # )
        # expected_cell_end = convert_unit(
            # np.array(
                # [[
                    # 5.9960617905,
                    # -0.0068118798,
                    # -0.0102043036,
                # ], [
                    # -0.0068116027,
                    # 6.0225574669,
                    # -0.0155044063,
                # ], [
                    # -0.0102048226,
                    # -0.0155046726,
                    # 6.0083072343,
                # ]]),
            # "angstrom"
        # )
        # self.assertEqual(simulation_cell.shape[0], 11)
        # self.assertTrue(np.array_equal(expected_cell_start, simulation_cell[0,:,:]))
        # self.assertTrue(np.array_equal(expected_cell_end, simulation_cell[-1,:,:]))


# #===============================================================================
# class TestElectronicStructureMethod(unittest.TestCase):

    # def test_mp2(self):
        # results = get_results("electronic_structure_method/mp2", "section_run")
        # result = results["electronic_structure_method"]
        # self.assertEqual(result, "MP2")

    # def test_dft_plus_u(self):
        # results = get_results("electronic_structure_method/dft_plus_u", "section_run")
        # result = results["electronic_structure_method"]
        # self.assertEqual(result, "DFT+U")

    # def test_rpa(self):
        # results = get_results("electronic_structure_method/rpa", "section_run")
        # result = results["electronic_structure_method"]
        # self.assertEqual(result, "RPA")


#===============================================================================
# class TestInputParser(unittest.TestCase):
    # """Tests that the parser can handle single-point calculations.
    # """

    # @classmethod
    # def setUpClass(cls):
        # cls.results = get_results("single_point", "section_run")
        # # cls.results.print_summary()

    # def test_x_cpmd_input_functional(self):
        # result = self.results["x_cpmd_input_DFT.FUNCTIONAL_options"]
        # self.assertEqual(result, "LDA")

    # def test_x_cpmd_input_cutoff(self):
        # result = self.results["x_cpmd_input_SYSTEM.CUTOFF_parameters"]
        # self.assertEqual(result, "70.0")

    # def test_x_cpmd_input_info(self):
        # result = self.results["x_cpmd_input_INFO_default_keyword"]
        # self.assertEqual(result, "isolated hydrogen molecule.\nsingle point calculation.")

    # def test_x_cpmd_input_atoms(self):
        # result = self.results["x_cpmd_input_ATOMS_default_keyword"]
        # self.assertEqual(result, "*H_MT_LDA.psp\nLMAX=S\n2\n4.371   4.000   4.000\n3.629   4.000   4.000")

    # def test_x_cpmd_input_optimize_wavefunction(self):
        # self.results["x_cpmd_section_input_CPMD.OPTIMIZE_WAVEFUNCTION"]


#===============================================================================
if __name__ == '__main__':
    suites = []
    suites.append(unittest.TestLoader().loadTestsFromTestCase(TestDFTEnergy))
    suites.append(unittest.TestLoader().loadTestsFromTestCase(TestDFTForce))
    suites.append(unittest.TestLoader().loadTestsFromTestCase(TestDFTGeoOpt))
    suites.append(unittest.TestLoader().loadTestsFromTestCase(TestXCFunctional))

    # suites.append(unittest.TestLoader().loadTestsFromTestCase(TestGeoOpt))
    # suites.append(unittest.TestLoader().loadTestsFromTestCase(TestInputParser))
    # suites.append(unittest.TestLoader().loadTestsFromTestCase(TestMD))
    # suites.append(unittest.TestLoader().loadTestsFromTestCase(TestMDTrajFormats))
    # suites.append(unittest.TestLoader().loadTestsFromTestCase(TestMDPrintSettings))
    # suites.append(unittest.TestLoader().loadTestsFromTestCase(TestPeriodicity))
    # suites.append(unittest.TestLoader().loadTestsFromTestCase(TestXCFunctional))
    alltests = unittest.TestSuite(suites)
    unittest.TextTestRunner(verbosity=0).run(alltests)
