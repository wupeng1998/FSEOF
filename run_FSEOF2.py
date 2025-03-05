#!/usr/bin/env python

from multiprocessing import freeze_support
from cobra.io import read_sbml_model
import cobra
import pandas as pd
import numpy as np
from cobra.flux_analysis import flux_variability_analysis
from scipy.optimize import curve_fit
import sys
import argparse
from FSEOF import FSEOF

#Parse arguments form the command line
parser = argparse.ArgumentParser(
    description="Identify gentetic targets for over-expression to increase the flux to a reaction of interest."
    )

parser.add_argument("sbmlFile", type=str, help="Path to the SBML model that should be scanned for over-expression targets. Needs to be a .xml file!")
parser.add_argument("biomassID", type=str, help="ID of the biomass reaction of the SBML file you provided")
parser.add_argument("reactionID", type=str, help="ID of the reaction that will be optimized by genetic engineering. Over-expression targets will be identified for this reaction of the SBML model")

parser.add_argument("--steps", type=int, action="store", default=30, help="Number of steps for the FSEOF algorithm. The default is 30. NOTE: This number should be decreased if flux variabillity is used.")
parser.add_argument("--useFVA", action="store_true", help="Changes the method for finding over-expression from flux balance analysis to flux variabillity analysis. This will significantly increase the runtime from a few minutes to several hours!")
parser.add_argument("--constrainBiomass", action="store_true", help="Constrains growth rate to 95%% of the theoretic maximum. This might improve the accuracy of the results, but can also lead to mathmatical infeasible solutions that will cause an error.")
parser.add_argument("--changeBiomassConstrain", type=float, action="store", default=0.95, help="If you would like to add an additional constrain to the growth rate, but not at 95%% of the theoretic maximum, use this option to specify at what percentage you want to set the constrain. NOTE: Percentages need to be passed as floats" )

args = parser.parse_args()

def parse_reaction(reaction):
    if pd.isna(reaction):
        return [], []
    
    # 先按反应符号分割
    if "<=>" in reaction:
        parts = reaction.split("<=>")
    elif "-->" in reaction:
        parts = reaction.split("-->")
    elif "<--" in reaction:
        parts = reaction.split("<--")
    else:
        return [], []
    
    # 进一步按 "+" 拆分
    left_compounds = [met.strip() for met in parts[0].split("+")]
    right_compounds = [met.strip() for met in parts[1].split("+")]
    
    # 按 "_" 拆分获取最终的代谢物名称
    left_elements = sorted(set([comp.split("_")[0] for comp in left_compounds]))
    right_elements = sorted(set([comp.split("_")[0] for comp in right_compounds]))
    
    return left_elements, right_elements

def main():
    
    f = FSEOF(args.sbmlFile, args.biomassID, args.reactionID)
    f.find_targets(args.steps, useFVA=args.useFVA, constrainBiomass=args.constrainBiomass, maxFluxCutoff=args.changeBiomassConstrain)
    f.addReactionData()

    # 定义辅助函数计算Biomass值
    def get_biomass_value(model, reaction_id):
        reaction = model.reactions.get_by_id(reaction_id)
        try:
            with model:
                reaction.bounds = (0,0)
                solution = model.optimize()
                if solution and solution.status == 'optimal':
                    return solution.objective_value
                else:
                    return 0.0
        except Exception as e:
            print(f"Error optimizing reaction {reaction_id}: {e}")
            return None

    filename = "AmplificationTargets_{reaction}.xlsx".format(reaction=args.reactionID)

    # 拆分 `q_slope`
    up_targets = f.targets[f.targets["q_slope"] >= 0].copy()
    down_targets = f.targets[f.targets["q_slope"] < 0].copy()

    # 重命名 `q_slope`
    up_targets = up_targets.rename(columns={"q_slope": "q_slope_up"})
    down_targets = down_targets.rename(columns={"q_slope": "q_slope_down"})

    # 删除 `Reaction_ID` 以 `EX_` 开头的行
    up_targets = up_targets[~up_targets["Reaction_ID"].str.startswith("EX_")]
    down_targets = down_targets[~down_targets["Reaction_ID"].str.startswith("EX_")]
    
    # 解析 Reaction 并删除左右代谢物相同的行
    up_targets[["Parsed_Left", "Parsed_Right"]] = up_targets["Reaction"].apply(lambda x: pd.Series(parse_reaction(x)))
    down_targets[["Parsed_Left", "Parsed_Right"]] = down_targets["Reaction"].apply(lambda x: pd.Series(parse_reaction(x)))
    
    up_targets = up_targets[up_targets["Parsed_Left"] != up_targets["Parsed_Right"]]
    down_targets = down_targets[down_targets["Parsed_Left"] != down_targets["Parsed_Right"]]
    
    # 删除辅助列
    up_targets = up_targets.drop(columns=["Parsed_Left", "Parsed_Right"])
    down_targets = down_targets.drop(columns=["Parsed_Left", "Parsed_Right"])

    # 添加Biomass列
    model = read_sbml_model(args.sbmlFile)
    # up_targets['Biomass'] = up_targets['Reaction_ID'].apply(lambda rid: get_biomass_value(f.model, rid))
    down_targets['Biomass'] = down_targets['Reaction_ID'].apply(lambda rid: get_biomass_value(model, rid))

    # 选择列（添加Biomass）
    if args.useFVA:
        columns = ["Reaction_ID", "q_slope_up", "l_sol", "q_slope_classifier", "l_sol_classifier", "reaction_class", "Reaction", "Compartments", "Genes"]
        columns_down = ["Reaction_ID", "q_slope_down", "l_sol", "q_slope_classifier", "l_sol_classifier", "reaction_class", "Reaction", "Compartments", "Biomass", "Genes"]
    else:
        columns = ["Reaction_ID", "q_slope_up", "Reaction", "Compartments", "Genes"]
        columns_down = ["Reaction_ID", "q_slope_down", "Reaction", "Compartments", "Biomass", "Genes"]

    up_targets = up_targets[columns]
    down_targets = down_targets[columns_down]

    # 按绝对值排序
    up_targets = up_targets.sort_values(by="q_slope_up", key=lambda x: x.abs(), ascending=False)
    down_targets = down_targets.sort_values(by="q_slope_down", key=lambda x: x.abs(), ascending=False)

    # 写入Excel
    with pd.ExcelWriter(filename) as writer:
        up_targets.to_excel(writer, sheet_name="up", index=False)
        down_targets.to_excel(writer, sheet_name="down", index=False)

    print(f"Results saved in {filename} with sheets: 'up' and 'down'")

if __name__ == "__main__":
    freeze_support()
    main()
