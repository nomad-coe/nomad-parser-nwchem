#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import pytest

from nomad.datamodel import EntryArchive
from nwchemparser.nwchem_parser import NWChemParser


@pytest.fixture(scope='module')
def parser():
    return NWChemParser()


def test_single_point(parser):
    archive = EntryArchive()
    parser.parse('tests/data/single_point.out', archive, None)

    sec_run = archive.section_run[0]
    assert sec_run.program_basis_set_type == 'gaussians'
    assert sec_run.program_version == '6.6'
    assert sec_run.x_nwchem_section_start_information[0].x_nwchem_ga_revision == '10594'

    sec_method = archive.section_run[0].section_method[0]
    assert sec_method.scf_max_iteration == 50
    assert sec_method.total_charge == 0
    assert len(sec_method.section_XC_functionals) == 2
    assert sec_method.section_XC_functionals[1].XC_functional_name == 'MGGA_C_TPSS'
    assert sec_method.section_XC_functionals[0].XC_functional_weight == 1.0

    assert archive.section_run[0].section_sampling_method[0].sampling_method == 'single_point'

    sec_scc = archive.section_run[0].section_single_configuration_calculation[0]
    assert sec_scc.energy_total.magnitude == pytest.approx(-3.332424186333889e-16)
    assert sec_scc.x_nwchem_energy_one_electron == pytest.approx(-5.35955587575652e-16)
    assert sec_scc.atom_forces[2][0].magnitude == pytest.approx(-4.9432341e-13)
    sec_scfs = sec_scc.section_scf_iteration
    assert len(sec_scfs) == 6
    assert sec_scfs[2].energy_total_scf_iteration.magnitude == pytest.approx(-3.33233301e-16)
    assert sec_scfs[5].time_scf_iteration.magnitude == 0.3
    assert sec_scfs[4].energy_change_scf_iteration.magnitude == pytest.approx(-7.45516347e-23)

    sec_system = archive.section_run[0].section_system[0]
    assert len(sec_system.atom_labels) == 3
    assert sec_system.atom_positions[0][2].magnitude == pytest.approx(-1.1817375e-11)


def test_geometry_optimization(parser):
    archive = EntryArchive()
    parser.parse('tests/data/geometry_optimization.out', archive, None)

    sec_methods = archive.section_run[0].section_method
    assert len(sec_methods) == 4
    assert sec_methods[1].scf_threshold_energy_change == pytest.approx(4.35974472e-24)

    sec_sccs = archive.section_run[0].section_single_configuration_calculation
    assert len(sec_sccs) == 4
    assert sec_sccs[0].energy_C.magnitude == pytest.approx(2.02746469e-16)
    assert sec_sccs[1].atom_forces[2][2].magnitude == pytest.approx(-2.20633015e-10)
    assert len(sec_sccs[2].section_scf_iteration) == 5

    sec_systems = archive.section_run[0].section_system
    assert sec_systems[0].atom_positions[1][2].magnitude == pytest.approx(5.6568542e-11)
    assert sec_systems[3].configuration_periodic_dimensions == [False, False, False]


def test_molecular_dynamics(parser):
    archive = EntryArchive()
    parser.parse('tests/data/molecular_dynamics.out', archive, None)

    sec_sccs = archive.section_run[0].section_single_configuration_calculation
    assert len(sec_sccs) == 6
    assert sec_sccs[2].energy_XC.magnitude == pytest.approx(-4.04565658e-17)
    assert sec_sccs[5].x_nwchem_section_qmd_step[0].x_nwchem_qmd_step_total_energy.magnitude == pytest.approx(-3.32745352e-16)
    assert sec_sccs[2].x_nwchem_section_qmd_step[0].x_nwchem_qmd_step_dipole[1] == pytest.approx(1.141435e-01)


def test_pw(parser):
    archive = EntryArchive()
    parser.parse('tests/data/pw.out', archive, None)

    sec_sccs = archive.section_run[0].section_single_configuration_calculation
    assert sec_sccs[1].energy_total.magnitude == pytest.approx(-8.89979631e-17)
    assert sec_sccs[1].spin_S2 == pytest.approx(2.0029484134705502)
    sec_scfs = sec_sccs[0].section_scf_iteration
    assert len(sec_scfs) == 5
    assert sec_scfs[3].energy_total_scf_iteration.magnitude == pytest.approx(-8.86651108e-17)