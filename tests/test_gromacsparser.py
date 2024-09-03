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
import numpy as np
import h5py
import json

from typing import (
    Any,
    Dict,
)

from nomad.datamodel import EntryArchive
from atomisticparsers.gromacs import GromacsParser
from simulationworkflowschema.molecular_dynamics import FreeEnergyCalculationParameters
from atomisticparsers.gromacs import GromacsLogParser


def approx(value, abs=0, rel=1e-6):
    return pytest.approx(value, abs=abs, rel=rel)


@pytest.fixture(scope='module')
def parser():
    return GromacsParser()


def test_md_verbose(parser):
    archive = EntryArchive()
    parser.parse('tests/data/gromacs/fe_test/md.log', archive, None)

    sec_run = archive.run[0]
    assert sec_run.program.version == '5.1.4'
    sec_control = sec_run.x_gromacs_section_control_parameters
    assert sec_control.x_gromacs_all_input_parameters.get('coulombtype') == 'pme'
    assert np.shape(
        sec_control.x_gromacs_all_input_parameters.get('qm-opts').get('deform')
    ) == (3, 3)

    sec_workflow = archive.workflow2
    assert sec_workflow.m_def.name == 'MolecularDynamics'
    sec_method = sec_workflow.method
    assert sec_method.thermodynamic_ensemble == 'NPT'
    assert sec_method.integrator_type == 'leap_frog'
    assert sec_method.integration_timestep.magnitude == 5e-16
    assert sec_method.integration_timestep.units == 'second'
    assert sec_method.n_steps == 20
    assert sec_method.coordinate_save_frequency == 20
    assert sec_method.thermodynamics_save_frequency == 5
    assert sec_method.thermostat_parameters[0].thermostat_type == 'berendsen'
    assert sec_method.thermostat_parameters[0].reference_temperature.magnitude == 298.0
    assert sec_method.thermostat_parameters[0].reference_temperature.units == 'kelvin'
    assert sec_method.thermostat_parameters[0].coupling_constant.magnitude == 5e-13
    assert sec_method.thermostat_parameters[0].coupling_constant.units == 'second'
    assert sec_method.barostat_parameters[0].barostat_type == 'berendsen'
    assert sec_method.barostat_parameters[0].coupling_type == 'isotropic'
    assert np.all(
        sec_method.barostat_parameters[0].reference_pressure.magnitude
        == [[100000.0, 0.0, 0.0], [0.0, 100000.0, 0.0], [0.0, 0.0, 100000.0]]
    )
    assert sec_method.barostat_parameters[0].reference_pressure.units == 'pascal'
    assert np.all(
        sec_method.barostat_parameters[0].coupling_constant.magnitude
        == [
            [1.0e-12, 1.0e-12, 1.0e-12],
            [1.0e-12, 1.0e-12, 1.0e-12],
            [1.0e-12, 1.0e-12, 1.0e-12],
        ]
    )
    assert sec_method.barostat_parameters[0].coupling_constant.units == 'second'
    assert np.all(
        sec_method.barostat_parameters[0].compressibility.magnitude
        == [
            [4.6e-10, 0.0e00, 0.0e00],
            [0.0e00, 4.6e-10, 0.0e00],
            [0.0e00, 0.0e00, 4.6e-10],
        ]
    )
    assert sec_method.barostat_parameters[0].compressibility.units == '1 / pascal'

    sec_sccs = sec_run.calculation
    assert len(sec_sccs) == 5
    assert sec_sccs[1].pressure_tensor[1][2].magnitude == approx(40267181.396484375)
    assert sec_sccs[3].pressure.magnitude == approx(-63926916.50390625)
    assert sec_sccs[3].temperature.magnitude == approx(291.80401611328125)
    assert sec_sccs[2].volume.magnitude == approx(1.505580043792725e-26)
    assert sec_sccs[2].density.magnitude == approx(1007.9478759765625)
    assert sec_sccs[2].enthalpy.magnitude == approx(-3.265048185525196e-17)
    assert sec_sccs[2].virial_tensor[2][2].magnitude == approx(1.1367756347656254e-19)
    assert len(sec_sccs[1].x_gromacs_thermodynamics_contributions) == 5
    assert sec_sccs[1].x_gromacs_thermodynamics_contributions[2].kind == '#Surf*SurfTen'
    assert sec_sccs[1].x_gromacs_thermodynamics_contributions[2].value == approx(
        2453.242431640625
    )
    assert len(sec_sccs[4].energy.x_gromacs_energy_contributions) == 12
    assert sec_sccs[-2].energy.x_gromacs_energy_contributions[1].kind == 'G96Angle'
    assert sec_sccs[-2].energy.x_gromacs_energy_contributions[
        1
    ].value.magnitude == approx(2.7314541413859703e-20)
    assert sec_sccs[0].energy.total.value.magnitude == approx(-3.2711273151685104e-17)
    assert sec_sccs[0].energy.electrostatic.value.magnitude == approx(
        -4.598738980683062e-17
    )
    assert sec_sccs[0].energy.electrostatic.short_range.magnitude == approx(
        -4.155359489335505e-17
    )
    assert sec_sccs[0].energy.electrostatic.long_range.magnitude == approx(
        -4.4337949134755685e-18
    )
    assert sec_sccs[-1].energy.van_der_waals.value.magnitude == approx(
        7.168028758853762e-18
    )
    assert sec_sccs[-1].energy.van_der_waals.short_range.magnitude == approx(
        7.37736630833841e-18
    )
    assert sec_sccs[-1].energy.van_der_waals.long_range.magnitude == approx(
        -1.2185287133832595e-19
    )
    assert sec_sccs[-1].energy.van_der_waals.correction.magnitude == approx(
        -8.74846781463227e-20
    )
    assert sec_sccs[0].energy.pressure_volume_work.value.magnitude == approx(
        1.5056965850293703e-21
    )

    assert sec_sccs[0].forces.total.value[5][2].magnitude == approx(
        -7.932968909721231e-10
    )

    sec_systems = sec_run.system
    assert len(sec_systems) == 2
    assert np.shape(sec_systems[0].atoms.positions) == (1516, 3)
    assert sec_systems[1].atoms.positions[800][1].magnitude == approx(2.4740036e-09)
    assert sec_systems[0].atoms.velocities[500][0].magnitude == approx(869.4773)
    assert sec_systems[1].atoms.lattice_vectors[2][2].magnitude == approx(2.469158e-09)
    assert sec_systems[0].atoms.bond_list[200, 0] == 289

    sec_method = sec_run.method
    assert len(sec_method) == 1
    assert len(sec_method[0].force_field.model[0].contributions) == 8
    assert sec_method[0].force_field.model[0].contributions[6].type == 'bond'
    assert sec_method[0].force_field.model[0].contributions[6].n_interactions == 1017
    assert sec_method[0].force_field.model[0].contributions[6].n_atoms == 2
    assert sec_method[0].force_field.model[0].contributions[6].atom_labels[10][0] == 'C'
    assert (
        sec_method[0].force_field.model[0].contributions[6].atom_indices[100, 1] == 141
    )
    assert sec_method[0].force_field.force_calculations.vdw_cutoff.magnitude == 1.2e-09
    assert sec_method[0].force_field.force_calculations.vdw_cutoff.units == 'meter'
    assert (
        sec_method[0].force_field.force_calculations.coulomb_type
        == 'particle_mesh_ewald'
    )
    assert sec_method[0].force_field.force_calculations.coulomb_cutoff.magnitude == 0.9
    assert sec_method[0].force_field.force_calculations.coulomb_cutoff.units == 'meter'
    assert (
        sec_method[
            0
        ].force_field.force_calculations.neighbor_searching.neighbor_update_frequency
        == 5
    )
    assert (
        sec_method[
            0
        ].force_field.force_calculations.neighbor_searching.neighbor_update_cutoff.magnitude
        == 9.000000000000001e-10
    )
    assert (
        sec_method[
            0
        ].force_field.force_calculations.neighbor_searching.neighbor_update_cutoff.units
        == 'meter'
    )


def test_md_edr(parser):
    archive = EntryArchive()
    parser.parse('tests/data/gromacs/fe_test/mdrun.out', archive, None)

    assert len(archive.run[0].calculation) == 5


def test_md_atomsgroup(parser):
    archive = EntryArchive()
    parser.parse(
        'tests/data/gromacs/polymer_melt/step4.0_minimization.log', archive, None
    )

    sec_run = archive.run[0]
    sec_systems = sec_run.system

    assert len(sec_systems[0].atoms_group) == 1
    assert len(sec_systems[0].atoms_group[0].atoms_group) == 100

    assert sec_systems[0].atoms_group[0].label == 'group_S1P1'
    assert sec_systems[0].atoms_group[0].type == 'molecule_group'
    assert sec_systems[0].atoms_group[0].index == 0
    assert sec_systems[0].atoms_group[0].composition_formula == 'S1P1(100)'
    assert sec_systems[0].atoms_group[0].n_atoms == 7200
    assert sec_systems[0].atoms_group[0].atom_indices[5] == 5
    assert sec_systems[0].atoms_group[0].is_molecule is False

    assert sec_systems[0].atoms_group[0].atoms_group[52].label == 'S1P1'
    assert sec_systems[0].atoms_group[0].atoms_group[52].type == 'molecule'
    assert sec_systems[0].atoms_group[0].atoms_group[52].index == 52
    assert (
        sec_systems[0].atoms_group[0].atoms_group[52].composition_formula == 'ETHOX(10)'
    )
    assert sec_systems[0].atoms_group[0].atoms_group[52].n_atoms == 72
    assert sec_systems[0].atoms_group[0].atoms_group[52].atom_indices[8] == 3752
    assert sec_systems[0].atoms_group[0].atoms_group[52].is_molecule is True

    assert (
        sec_systems[0].atoms_group[0].atoms_group[76].atoms_group[0].label
        == 'group_ETHOX'
    )
    assert (
        sec_systems[0].atoms_group[0].atoms_group[76].atoms_group[0].type
        == 'monomer_group'
    )
    assert sec_systems[0].atoms_group[0].atoms_group[76].atoms_group[0].index == 0
    assert (
        sec_systems[0].atoms_group[0].atoms_group[76].atoms_group[0].composition_formula
        == 'ETHOX(10)'
    )
    assert sec_systems[0].atoms_group[0].atoms_group[76].atoms_group[0].n_atoms == 72
    assert (
        sec_systems[0].atoms_group[0].atoms_group[76].atoms_group[0].atom_indices[5]
        == 5477
    )
    assert (
        sec_systems[0].atoms_group[0].atoms_group[76].atoms_group[0].is_molecule
        is False
    )

    assert (
        sec_systems[0]
        .atoms_group[0]
        .atoms_group[76]
        .atoms_group[0]
        .atoms_group[7]
        .label
        == 'ETHOX'
    )
    assert (
        sec_systems[0].atoms_group[0].atoms_group[76].atoms_group[0].atoms_group[7].type
        == 'monomer'
    )
    assert (
        sec_systems[0]
        .atoms_group[0]
        .atoms_group[76]
        .atoms_group[0]
        .atoms_group[7]
        .index
        == 7
    )
    assert (
        sec_systems[0]
        .atoms_group[0]
        .atoms_group[76]
        .atoms_group[0]
        .atoms_group[7]
        .composition_formula
        == 'C(2)H(4)O(1)'
    )
    assert (
        sec_systems[0]
        .atoms_group[0]
        .atoms_group[76]
        .atoms_group[0]
        .atoms_group[7]
        .n_atoms
        == 7
    )
    assert (
        sec_systems[0]
        .atoms_group[0]
        .atoms_group[76]
        .atoms_group[0]
        .atoms_group[7]
        .atom_indices[5]
        == 5527
    )
    assert (
        sec_systems[0]
        .atoms_group[0]
        .atoms_group[76]
        .atoms_group[0]
        .atoms_group[7]
        .is_molecule
        is False
    )


def test_geometry_optimization(parser):
    archive = EntryArchive()
    parser.parse(
        'tests/data/gromacs/polymer_melt/step4.0_minimization.log', archive, None
    )

    sec_workflow = archive.workflow2

    assert sec_workflow.method.type == 'atomic'
    assert sec_workflow.method.method == 'steepest_descent'
    assert sec_workflow.method.convergence_tolerance_force_maximum.magnitude == approx(
        6.02214076e38
    )
    assert sec_workflow.method.convergence_tolerance_force_maximum.units == 'newton'
    assert sec_workflow.results.final_force_maximum.magnitude == approx(
        1.303670442204273e38
    )
    assert sec_workflow.results.final_force_maximum.units == 'newton'
    assert sec_workflow.results.optimization_steps == 12
    assert sec_workflow.method.optimization_steps_maximum == 5000
    assert len(sec_workflow.results.energies) == 11
    assert sec_workflow.results.energies[2].magnitude == approx(8.244726173423347e-17)
    assert sec_workflow.results.energies[2].units == 'joule'
    assert len(sec_workflow.results.steps) == 11
    assert sec_workflow.results.steps[4] == 5000


def test_integrator_sd(parser):
    archive = EntryArchive()
    parser.parse(
        'tests/data/gromacs/water_AA_ENUM_tests/integrator-sd/md.log', archive, None
    )

    sec_run = archive.run[0]
    # assert sec_run.program.version == "2018.6"

    sec_workflow = archive.workflow2
    assert sec_workflow.m_def.name == 'MolecularDynamics'
    sec_method = sec_workflow.method
    assert sec_method.thermodynamic_ensemble == 'NVT'
    assert sec_method.integrator_type == 'langevin_goga'
    assert sec_method.thermostat_parameters[0].thermostat_type == 'langevin_goga'
    assert sec_method.thermostat_parameters[0].reference_temperature.magnitude == 298.0
    assert sec_method.thermostat_parameters[0].coupling_constant.magnitude == 5e-13


def test_integrator_mdvv(parser):
    archive = EntryArchive()
    parser.parse(
        'tests/data/gromacs/water_AA_ENUM_tests/integrator-mdvv/md.log', archive, None
    )

    sec_run = archive.run[0]
    # assert sec_run.program.version == "2018.6"

    sec_workflow = archive.workflow2
    assert sec_workflow.m_def.name == 'MolecularDynamics'
    sec_method = sec_workflow.method
    assert sec_method.thermodynamic_ensemble == 'NVE'
    assert sec_method.integrator_type == 'velocity_verlet'


def test_integrator_bd(parser):
    archive = EntryArchive()
    parser.parse(
        'tests/data/gromacs/water_AA_ENUM_tests/integrator-bd/md.log', archive, None
    )

    sec_run = archive.run[0]
    # assert sec_run.program.version == "2018.6"

    sec_workflow = archive.workflow2
    assert sec_workflow.m_def.name == 'MolecularDynamics'
    sec_method = sec_workflow.method
    assert sec_method.thermodynamic_ensemble == 'NVE'
    assert sec_method.integrator_type == 'brownian'


# TODO test for andersen thermostat? It's not clear how to run this at the moment or if it is deprecated in newer versions of Gromacs.


def test_integrator_md_thermostat_vrescale(parser):
    archive = EntryArchive()
    parser.parse(
        'tests/data/gromacs/water_AA_ENUM_tests/integrator-md/thermostat-vrescale/md.log',
        archive,
        None,
    )

    sec_run = archive.run[0]
    assert sec_run.program.version == '2018.6'

    sec_workflow = archive.workflow2
    assert sec_workflow.m_def.name == 'MolecularDynamics'
    sec_method = sec_workflow.method
    assert sec_method.thermodynamic_ensemble == 'NVT'
    assert sec_method.integrator_type == 'leap_frog'
    assert sec_method.thermostat_parameters[0].thermostat_type == 'velocity_rescaling'
    assert sec_method.thermostat_parameters[0].reference_temperature.magnitude == 298.0
    assert sec_method.thermostat_parameters[0].coupling_constant.magnitude == 5e-13


def test_integrator_md_thermostat_nosehoover_barostat_parrinellorahman(parser):
    archive = EntryArchive()
    parser.parse(
        'tests/data/gromacs/water_AA_ENUM_tests/integrator-md/thermostat-nosehoover_barostat-parrinellorahman/md.log',
        archive,
        None,
    )

    sec_run = archive.run[0]
    assert sec_run.program.version == '2018.6'

    sec_workflow = archive.workflow2
    assert sec_workflow.m_def.name == 'MolecularDynamics'
    sec_method = sec_workflow.method
    assert sec_method.thermodynamic_ensemble == 'NPT'
    assert sec_method.integrator_type == 'leap_frog'
    assert sec_method.thermostat_parameters[0].thermostat_type == 'nose_hoover'
    assert sec_method.thermostat_parameters[0].reference_temperature.magnitude == 298.0
    assert sec_method.thermostat_parameters[0].coupling_constant.magnitude == 5e-13
    assert sec_method.barostat_parameters[0].barostat_type == 'parrinello_rahman'
    assert sec_method.barostat_parameters[0].coupling_type == 'isotropic'
    assert np.all(
        sec_method.barostat_parameters[0].reference_pressure.magnitude
        == [[100000.0, 0.0, 0.0], [0.0, 100000.0, 0.0], [0.0, 0.0, 100000.0]]
    )
    assert np.all(
        sec_method.barostat_parameters[0].coupling_constant.magnitude
        == [
            [5.0e-12, 5.0e-12, 5.0e-12],
            [5.0e-12, 5.0e-12, 5.0e-12],
            [5.0e-12, 5.0e-12, 5.0e-12],
        ]
    )
    assert np.all(
        sec_method.barostat_parameters[0].compressibility.magnitude
        == [
            [7.4e-10, 0.0e00, 0.0e00],
            [0.0e00, 7.4e-10, 0.0e00],
            [0.0e00, 0.0e00, 7.4e-10],
        ]
    )


def test_free_energy_calculations(parser):
    archive = EntryArchive()
    parser.parse(
        'tests/data/gromacs/free_energy_calculations/alchemical_transformation_single_run/fep_run-7.log',
        archive,
        None,
    )

    import h5py

    def get_dataset(filename_with_path):
        try:
            # Split the filename and dataset path
            filename, dataset_path = filename_with_path.split('#', 1)

            # Open the HDF5 file in read mode
            with h5py.File(filename, 'r') as file:
                # Access the dataset using the provided path
                dataset = file[dataset_path]
                data = dataset[()]

            return data

        except (ValueError, KeyError) as e:
            # Handle potential errors (e.g., invalid input or dataset not found)
            print(f'Error: {e}')
            return None

    sec_workflow = archive.workflow2
    sec_method = sec_workflow.method.free_energy_calculation_parameters[0]
    sec_results = sec_workflow.results.free_energy_calculations[0]

    assert sec_method.type == 'alchemical'
    sec_lambdas = sec_method.lambdas
    assert len(sec_lambdas) == 7
    assert sec_lambdas[2].type == 'vdw'
    assert sec_lambdas[2].value[2] == 0.2
    assert sec_lambdas[-1].type == 'temperature'
    assert sec_lambdas[-1].value[2] == 0.0
    assert sec_method.lambda_index == 7
    assert sec_method.atom_indices.shape == (1,)
    assert sec_method.atom_indices[0] == 0
    assert sec_method.initial_state_vdw is True
    assert sec_method.final_state_vdw is False
    assert sec_method.initial_state_coloumb is False
    assert sec_method.final_state_coloumb is False
    assert sec_method.initial_state_bonded is True
    assert sec_method.final_state_bonded is True

    assert sec_results.n_frames == 5001
    assert sec_results.n_states == 11
    assert sec_results.lambda_index == 7
    assert len(sec_results.times) == 5001
    assert sec_results.times.to('ps')[10].magnitude == approx(2.0)
    assert sec_results.value_unit == 'kilojoule'
    # assert isinstance(sec_results.method_ref, FreeEnergyCalculationParameters)
    # TODO add testing of hdf5 references in sec_results ('value_total_energy_magnitude', 'value_total_energy_derivative_magnitude', 'value_total_energy_differences_magnitude', 'value_PV_energy_magnitude') to NOMAD testing


@pytest.mark.parametrize(
    'path, input_log_fnm, result_json_fnm',
    [
        ('tests/data/gromacs/input_parameters', 'test.log', 'input_parameters.json'),
        (
            'tests/data/gromacs/free_energy_calculations/alchemical_transformation_single_run',
            'fep_run-7.log',
            'input_parameters.json',
        ),
    ],
)
def test_str_to_input_parameters(path: str, input_log_fnm: str, result_json_fnm: str):
    """_summary_

    Args:
        input (str): _description_
        result (Dict[Any]): _description_
    """

    def assert_dict_equal(d1, d2):
        """
        Recursively assert that two dictionaries are equal.

        Args:
            d1 (dict): First dictionary to compare.
            d2 (dict): Second dictionary to compare.
        """

        assert isinstance(d1, dict), f'Expected dict, got {type(d1)}'
        assert isinstance(d2, dict), f'Expected dict, got {type(d2)}'
        assert d1.keys() == d2.keys(), f'Keys mismatch: {d1.keys()} != {d2.keys()}'

        for key in d1:
            if isinstance(d1[key], dict) and isinstance(d2[key], dict):
                assert_dict_equal(d1[key], d2[key])
            else:
                if isinstance(d1[key], (str, bool)):
                    assert (
                        d1[key] == d2[key]
                    ), f"Value mismatch for key '{key}': {d1[key]} != {d2[key]}"
                elif isinstance(d1[key], np.ndarray):
                    assert np.isclose(
                        d1[key], d2[key]
                    ).all(), f"Value mismatch for key '{key}': {d1[key]} != {d2[key]}"
                elif abs(d1[key]) == float('inf'):
                    assert 'inf' == d2[key] if d1[key] > 0 else '-inf' == d2[key]
                else:
                    assert d1[key] == approx(
                        d2[key]
                    ), f"Value mismatch for key '{key}': {d1[key]} != {d2[key]}"

    log_parser = GromacsLogParser()
    log_parser.mainfile = f'{path}/{input_log_fnm}'
    parsed_parameters = log_parser.get('input_parameters')
    parsed_parameters = {
        key.replace('_', '-'): val.lower() if isinstance(val, str) else val
        for key, val in parsed_parameters.items()
    }

    with open(f'{path}/{result_json_fnm}', 'r', encoding='utf-8') as f:
        result = json.load(f)
    __ = result.pop('mdp_unique_params')

    assert_dict_equal(parsed_parameters, result)
