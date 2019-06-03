# -*- coding: utf-8 -*-
from __future__ import absolute_import
from six.moves import range
__version__ = "0.3.1"

import re
from math import isnan
from aiida.parsers.parser import Parser
from aiida.orm.data.parameter import ParameterData
from aiida_raspa.calculations import RaspaCalculation
from aiida.parsers.exceptions import OutputParsingError

import pandas as pd
from pathlib import Path

float_base = float


def float(number):  # pylint: disable=redefined-builtin
    number = float_base(number)
    return number if not isnan(number) else None


KELVIN_TO_KJ_PER_MOL = float(8.314464919 / 1000.0)  #exactly the same as Raspa

block1_list = [
    (re.compile("Average Volume:"), "cell_volume", (1, 2, 4)),
    (re.compile("Average Pressure:"), "pressure", (1, 2, 4)),
    (re.compile("Average temperature:"), "temperature", (1, 2, 4)),
    (re.compile("Average Density:"), "adsorbate_density", (1, 2, 4)),
    (re.compile("Total energy:$"), "total_energy", (1, 2, 4)),
    (re.compile("Average Heat Capacity"), "framework_heat_capacity", (1, 2,
                                                                      4)),
    (re.compile("Heat of desorption:"), "heat_of_desorption", (1, 4, 3)),
    (re.compile("Enthalpy of adsorption:"), "enthalpy_of_adsorption", (1, 4,
                                                                       3)),
    (re.compile("Tail-correction energy:$"), "tail_correction_energy", (1, 2,
                                                                        4)),
]


# pylint: disable=too-many-arguments
def parse_block1(fl, result_dict, prop, value=1, units=2, dev=4):
    """Parse block that looks as follows:
    Average Volume:
    =================
        Block[ 0]        12025.61229 [A^3]
        Block[ 1]        12025.61229 [A^3]
        Block[ 2]        12025.61229 [A^3]
        Block[ 3]        12025.61229 [A^3]
        Block[ 4]        12025.61229 [A^3]
        ------------------------------------------------------------------------------
        Average          12025.61229 [A^3] +/-            0.00000 [A^3]
    """
    for line in fl:
        if 'Average' in line:
            result_dict[prop + '_average'] = float(line.split()[value])
            result_dict[prop + '_units'] = line.split()[units].translate(
                None, '[](){}')
            result_dict[prop + '_dev'] = float(line.split()[dev])
            return


energy_block_list = [
    (re.compile("Average Adsorbate-Adsorbate energy:"), 'ads_ads'),
    (re.compile("Average Host-Adsorbate energy:"), 'host_ads'),
]


def parse_block_energy(fl, res_dict, prop):
    """Parse block that looks as follows:
    Average Adsorbate-Adsorbate energy:
    ===================================
        Block[ 0] -443.23204         Van der Waals: -443.23204         Coulomb: 0.00000            [K]
        Block[ 1] -588.20205         Van der Waals: -588.20205         Coulomb: 0.00000            [K]
        Block[ 2] -538.43355         Van der Waals: -538.43355         Coulomb: 0.00000            [K]
        Block[ 3] -530.00960         Van der Waals: -530.00960         Coulomb: 0.00000            [K]
        Block[ 4] -484.15106         Van der Waals: -484.15106         Coulomb: 0.00000            [K]
        ------------------------------------------------------------------------------
        Average   -516.80566         Van der Waals: -516.805659        Coulomb: 0.00000            [K]
              +/- 98.86943                      +/- 98.869430               +/- 0.00000            [K]
    """
    for line in fl:
        if 'Average' in line:
            res_dict[prop + '_total_energy_unit'] = 'kJ/mol'
            res_dict[prop + '_vdw_energy_unit'] = 'kJ/mol'
            res_dict[prop + '_coulomb_energy_unit'] = 'kJ/mol'
            res_dict[prop + '_total_energy_average'] = float(
                line.split()[1]) * KELVIN_TO_KJ_PER_MOL
            res_dict[prop + '_vdw_energy_average'] = float(
                line.split()[5]) * KELVIN_TO_KJ_PER_MOL
            res_dict[prop + '_coulomb_energy_average'] = float(
                line.split()[7]) * KELVIN_TO_KJ_PER_MOL
        if '+/-' in line:
            res_dict[prop + '_total_energy_dev'] = float(
                line.split()[1]) * KELVIN_TO_KJ_PER_MOL
            res_dict[prop + '_vdw_energy_dev'] = float(
                line.split()[3]) * KELVIN_TO_KJ_PER_MOL
            res_dict[prop + '_coulomb_energy_dev'] = float(
                line.split()[5]) * KELVIN_TO_KJ_PER_MOL
            return


lines_with_component_list = [
    (re.compile(" Average chemical potential: "), "chemical_potential"),
    (re.compile(" Average Henry coefficient: "), "henry_coefficient"),
    (re.compile(" Average  <U_gh>_1-<U_h>_0:"), "adsorption_energy_widom"),
    (re.compile(" Average Widom Rosenbluth-weight:"),
     "widom_rosenbluth_factor"),
]


def parse_lines_with_component(res_components, components, line, prop):
    # self.logger.info("analysing line: {}".format(line))
    for i, c in enumerate(components):
        if '[' + c + ']' in line:
            words = line.split()
            res_components[i][prop + '_units'] = words[-1].translate(
                None, '[](){}')
            res_components[i][prop + '_dev'] = float(words[-2])
            res_components[i][prop + '_average'] = float(words[-4])


def parse_rdfs(rdf_file_list):
    rdf_dict = {}
    for file in rdf_file_list:
        name = Path(file).stem.strip('RDF_')
        df = pd.read_csv(file, comment='#', header=None, sep='\s+')
        r = df.iloc[:, 1].values
        rdf = df.iloc[:, 2].values
        rdf_dict[name] = {'r': r, 'rdf': rdf}
    return rdf_dict


def parse_performance_line(line):
    """
    :param line: string in format
        Component [tip4pew] total tried: 386.000000 succesfull growth: 3.000000 (0.777202 [%]) accepted: 3.000000 (0.777202 [%])
        that will be parsed
    :return:
    """
    try:
        parts = line.split()
        total_tried = parts[4]
        successfull_growth_total = parts[7]
        successfull_growth_ratio = parts[8].strip('(')
        accepted_total = parts[-3]
        accepted_ratio = parts[-2].strip('(')
        return {
            'total_tried': int(float(total_tried)),
            'successfull_growth_total': int(float(successfull_growth_total)),
            'successfull_growth_ratio': float(successfull_growth_ratio),
            'accepted_total': int(float(accepted_total)),
            'accepted_ratio': float(accepted_ratio)
        }
    except Exception:
        return {}


def parse_performance_block(lines):
    """
    :param lines: list of strings in format
        	total        1.000000 1.000000 2.000000
            succesfull   1.000000 1.000000 1.000000
            accepted   1.000000 1.000000 0.500000
            displacement 0.487500 0.487500 0.325000
    :return:
    """
    import numpy as np
    try:
        totals = [int(float(i)) for i in lines[4].split()[1:]]
        successfull = [int(float(i)) for i in lines[3].split()[1:]]
        acceptance_ratio = [float(i) for i in lines[2].split()[1:]]
        drift = [float(i) for i in lines[1].split()[1:]]
        return {
            'total': totals,
            'successfull': successfull,
            'acceptance_ratio': acceptance_ratio,
            'drift': drift,
            'acceptance_ratio_mean': np.mean(acceptance_ratio)
        }
    except Exception:
        return {}


def parse_performance_mc(f):
    """
    Parse for some performance metrics of the MC moves
    :param f: file as lines
    :return: dictionary with efficiency of MC efficiencies
    """
    efficiencies_dict = {}
    # read from end for efficiency:
    for i, line in enumerate(f[::-1]):
        if ('Performance of the Reinsertion move:' in line) and ('OFF' not in line):
            efficiencies_dict['reinsertion'] = parse_performance_line(
                f[::-1][i - 2])
        if ('Performance of the swap deletion move:' in line) and ('OFF' not in line):
            efficiencies_dict['deletion'] = parse_performance_line(f[::-1][i -
                                                                           2])
        if ('Performance of the swap addition move:' in line) and ('OFF' not in line):
            efficiencies_dict['addition'] = parse_performance_line(f[::-1][i -
                                                                           2])
        if ('Performance of the rotation move:' in line)  and ('OFF' not in line):
            efficiencies_dict['rotation'] = parse_performance_block(
                f[::-1][i - 7:i - 2])

        if ('Performance of the translation move:' in line) and ('OFF' not in line):
            efficiencies_dict['translation'] = parse_performance_block(
                f[::-1][i - 7:i - 2])

        if 'Monte-Carlo moves statistics' in line:
            break

    return efficiencies_dict


class RaspaParser(Parser):
    """Parser for the output of RASPA."""

    # --------------------------------------------------------------------------
    def __init__(self, calc):
        """Initialize the instance of RaspaParser."""
        super(RaspaParser, self).__init__(calc)

        # check for valid input
        if not isinstance(calc, RaspaCalculation):
            raise OutputParsingError("Input calc must be a RaspaCalculation")

    # --------------------------------------------------------------------------
    def parse_with_retrieved(self, retrieved):
        """Receives in input a dictionary of retrieved nodes. Does all the logic here."""
        out_folder = retrieved['retrieved']
        new_nodes_list = []
        self._parse_stdout(out_folder, new_nodes_list)
        return True, new_nodes_list

    # pylint: disable=too-many-locals, too-many-arguments, too-many-statements, too-many-branches
    def _parse_stdout(self, out_folder, new_nodes_list):
        fn = None
        fs = out_folder.get_folder_list()
        for f in fs:
            if f.endswith('.data'):
                fn = f
        if fn is None:
            raise OutputParsingError(
                "Calculation did not produce an output file. Please make sure that it run "
                "correctly")

        res_per_component = []
        component_names = []
        inp_params = self._calc.inp.parameters.get_dict()
        ncomponents = len(inp_params['Component'])
        for i in range(ncomponents):
            res_per_component.append({})
            component_names.append(inp_params['Component'][i]['MoleculeName'])
        # self.logger.info("list of components: {}".format(component_names))

        output_abs_path = out_folder.get_abs_path(fn)
        result_dict = {'exceeded_walltime': False}
        framework_density = re.compile("Framework Density:")
        num_of_molec = re.compile("Number of molecules:$")

        rdfs = [out_folder.get_abs_path(f) for f in fs if f.startswith('RDF')]
        if len(rdfs) > 0:
            rdf_dict = parse_rdfs(rdfs)
        else:
            rdf_dict = None

        result_dict['rdfs'] = rdf_dict

        with open(output_abs_path, "r") as f:
            # 1st parsing part
            icomponent = 0
            res_cmp = res_per_component[0]
            warnings = []
            for line in f:
                # TODO maybe change for parse_line?
                if "WARNING" in line:
                    self.logger.warning(line)
                    warnings.append(line.split('WARNING:')[-1])
                if "Conversion factor molecules/unit cell -> mol/kg:" in line:
                    res_cmp['conversion_factor_molec_uc_to_mol_kg'] = float(
                        line.split()[6])
                    res_cmp[
                        'conversion_factor_molec_uc_to_mol_kg_unit'] = "(mol/kg)/(molec/uc)"
                if "Conversion factor molecules/unit cell -> gr/gr:" in line:
                    res_cmp['conversion_factor_molec_uc_to_gr_gr'] = float(
                        line.split()[6])
                    res_cmp[
                        'conversion_factor_molec_uc_to_gr_gr_unit'] = "(gr/gr)/(molec/uc)"
                if "Conversion factor molecules/unit cell -> cm^3 STP/gr:" in line:
                    res_cmp['conversion_factor_molec_uc_to_cm3stp_gr'] = float(
                        line.split()[7])
                    res_cmp[
                        'conversion_factor_molec_uc_to_cm3stp_gr_unit'] = "(cm^3_STP/gr)/(molec/uc)"
                if "Conversion factor molecules/unit cell -> cm^3 STP/cm^3:" in line:
                    res_cmp[
                        'conversion_factor_molec_uc_to_cm3stp_cm3'] = float(
                            line.split()[7])
                    res_cmp[
                        'conversion_factor_molec_uc_to_cm3stp_cm3_unit'] = "(cm^3_STP/cm^3)/(molec/uc)"
                if "MolFraction:" in line:
                    res_cmp['mol_fraction'] = float(line.split()[1])
                    res_cmp['mol_fraction_unit'] = "-"
                if "Partial pressure:" in line:
                    res_cmp['partial_pressure'] = float(line.split()[2])
                    res_cmp['partial_pressure_unit'] = "Pa"
                if "Partial fugacity:" in line:
                    res_cmp['partial_fugacity'] = float(line.split()[2])
                    res_cmp['partial_fugacity_unit'] = "Pa"
                    icomponent += 1
                    if icomponent < ncomponents:
                        res_cmp = res_per_component[icomponent]
                    else:
                        break
            # end of the 1st parsing part

            # 2nd parsing part
            for line in f:
                for parse in block1_list:
                    if parse[0].match(line):
                        parse_block1(f, result_dict, parse[1], *parse[2])
                        for i, cmpnt in enumerate(component_names):
                            # I assume here that properties per component are present right in the next line. The order of
                            # properties per molecule is the same as the order of molecules in the input file. So if component
                            # name was not found in the next line, I break the loop immidiately as there is no reason to
                            # continue it
                            line = next(f)
                            if cmpnt in line:
                                parse_block1(f, res_per_component[i], parse[1],
                                             *parse[2])
                            else:
                                break

                        continue  # no need to perform further checks, property has been found already
                for parse in energy_block_list:
                    if parse[0].match(line):
                        parse_block_energy(f, result_dict, prop=parse[1])
                        continue  # no need to perform further checks, property has been found already
                if framework_density.match(line) is not None:
                    result_dict['framework_density'] = line.split()[2]
                    result_dict['framework_density_units'] = line.split(
                    )[3].translate(None, '[](){}')

                elif num_of_molec.match(line) is not None:
                    break  # this stops the cycle
            # end of the 2nd parsing part

            # 3rd parsing part
            icomponent = 0
            for line in f:
                # TODO: change for parse_line?
                if 'Average loading absolute [molecules/unit cell]' in line:
                    res_per_component[icomponent][
                        'loading_absolute_average'] = float(line.split()[5])
                    res_per_component[icomponent][
                        'loading_absolute_dev'] = float(line.split()[7])
                    res_per_component[icomponent][
                        'loading_absolute_units'] = 'molecules/unit cell'
                elif 'Average loading excess [molecules/unit cell]' in line:
                    res_per_component[icomponent][
                        'loading_excess_average'] = float(line.split()[5])
                    res_per_component[icomponent][
                        'loading_excess_dev'] = float(line.split()[7])
                    res_per_component[icomponent][
                        'loading_excess_units'] = 'molecules/unit cell'
                    icomponent += 1
                if icomponent >= ncomponents:
                    break
            # end of the 3rd parsing part

            # 4th parsing part
            for line in f:
                for to_parse in lines_with_component_list:
                    if to_parse[0].search(line):
                        parse_lines_with_component(res_per_component,
                                                   component_names, line,
                                                   to_parse[1])
            # end of the 4th parsing part

            result_dict['warnings'] = warnings

        with open(output_abs_path, "r") as f:
            mc_performance = parse_performance_mc(f.readlines())
            result_dict['mc_move_statistics'] = mc_performance

        pair = (self.get_linkname_outparams(), ParameterData(dict=result_dict))
        new_nodes_list.append(pair)
        for i, item in enumerate(res_per_component):
            pair = ('component_' + str(i), ParameterData(dict=item))
            new_nodes_list.append(pair)


# EOF
