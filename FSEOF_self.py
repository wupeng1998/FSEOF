import cobra
from cobra.flux_analysis import pfba
import pandas as pd
import argparse
from cobra.io import read_sbml_model

def fseof_analysis(model_path, biomass_reaction, objective_reaction):
    """
    进行FSEOF分析并输出强化、弱化和敲除靶点到Excel文件。
    
    参数:
        model_path (str): SBML模型文件路径
        objective_reaction (str): 目标反应ID
        biomass_reaction (str): 生物量反应ID
    """
    # 读取模型
    model = read_sbml_model(model_path)
    
    # 设置目标反应和生物量反应
    model.objective = objective_reaction
    biomass_reaction_obj = model.reactions.get_by_id(biomass_reaction)
    
    # 计算最大理论通量
    pfba_solution = pfba(model)
    max_theoretical_flux = pfba_solution.fluxes[objective_reaction]
    print(f"最大理论通量: {max_theoretical_flux}")
    
    # 设置通量水平
    max_enforced_flux = 0.9
    max_flux = max_theoretical_flux * max_enforced_flux
    number_of_results = 10
    levels = [ (i + 1) * (max_flux) / number_of_results for i in range(number_of_results)]
    print(f"通量水平: {levels}")
    
    # 计算不同通量水平下的代谢网络，利用with上下文隔离模型修改
    df_list = []
    for f in levels:
        with model:
            model.reactions.get_by_id(objective_reaction).bounds = (f, f)
            model.objective = biomass_reaction_obj
            sol = pfba(model)
            # 转换为DataFrame并去除'reduced_costs'列
            sol_df = sol.to_frame()
            sol_df.drop('reduced_costs', axis=1, inplace=True)
            df_list.append(sol_df)
    
    # 合并不同通量水平下的结果，列名为1~10
    df_new = pd.concat(df_list, axis=1)
    df_new.columns = list(range(1, number_of_results+1))
    df_new.fillna(value=0, inplace=True)
    
    # --- 优化数据处理部分 ---
    # 1. 去除所有通量绝对值均<=1e-5的行（向量化操作）
    mask_all_zero = (df_new.abs() <= 1e-5).all(axis=1)
    df_new = df_new.loc[~mask_all_zero]
    
    # 2. 删除指定的反应行
    remove_rea = ['GLCtex_copy1', 'NH4tex', 'NH4tpp', 'H2Otex', 'H2Otpp',
                  'EX_h2o_e', 'O2tpp', 'O2tex', 'EX_o2_e', 'CO2tpp', 'CO2tex',
                  'Htex', 'EX_h_e', 'ADD_H_c-tex', 'ADD_EX_h_c']
    df_new.drop(index=remove_rea, errors='ignore', inplace=True)
    
    # 3. 批量调整权重
    weight = ['ATPM', 'TPI', 'ENO', 'PGM', 'PGK', 'GLCtex_copy1', 'NH4tex', 'NH4tpp',
              'EX_nh4_e', 'H2Otex', 'H2Otpp', 'EX_h2o_e', 'O2tpp', 'O2tex',
              'EX_o2_e', 'CO2tpp', 'CO2tex', 'EX_co2_e', 'Htex', 'EX_h_e',
              'ADD_H_c-tex', 'ADD_EX_h_c']
    lose_weight = ['F6PA', 'FBA3', 'EDA', 'XYLI2']
    flux_cols = list(range(1, number_of_results+1))
    
    # 注意：df_new的索引为反应ID（字符串）
    df_new.loc[df_new.index.isin(weight), flux_cols] = df_new.loc[df_new.index.isin(weight), flux_cols] * 1000
    df_new.loc[df_new.index.isin(lose_weight), flux_cols] = df_new.loc[df_new.index.isin(lose_weight), flux_cols] / 50
    
    df_new = df_new.round(5)
    
    # 插入反应方程式列（利用列表解析，速度较快）
    reaction_equations = []
    for rxn_id in df_new.index:
        try:
            reaction_equation = model.reactions.get_by_id(rxn_id).reaction
        except Exception:
            reaction_equation = ""
        reaction_equations.append(reaction_equation)
    df_new.insert(0, 'reactions', reaction_equations)
    
    # 计算额外指标：初始值、末值、最大、最小绝对值（向量化计算）
    df_new['abs_initial'] = df_new[1].abs()
    df_new['abs_end'] = df_new[number_of_results].abs()
    df_new['max_abs'] = df_new[flux_cols].abs().max(axis=1)
    df_new['min_abs'] = df_new[flux_cols].abs().min(axis=1)
    
    # 筛选强化靶点（upregulated）
    # 条件：初末符号一致且末端绝对值增大，同时末端接近最大值
    mask_up = (df_new[1] * df_new[number_of_results] >= 0) & \
              (df_new['abs_end'] > df_new['abs_initial']) & (
                  (df_new['max_abs'] == df_new['abs_end']) | 
                  (((df_new['max_abs'] - df_new['abs_end']) / df_new['max_abs']) < 0.1)
              )
    upregulated = df_new.index[mask_up]
    
    # 筛选弱化靶点（downregulated）
    mask_down = (df_new[1] * df_new[number_of_results] >= 0) & \
                (df_new['abs_end'] < df_new['abs_initial']) & \
                (df_new['min_abs'] == df_new['abs_end'])
    downregulated = df_new.index[mask_down]
    
    # 筛选敲除靶点（knockout）：末端通量为0且初始非0
    mask_ko = (df_new[number_of_results] == 0) & (df_new[1] != 0)
    knockout = df_new.index[mask_ko]
    
    # 辅助函数：构造输出数据（对每个反应调用model.reactions.get_by_id获取反应详细信息）
    def build_output_data(rxn_ids):
        data = []
        for rxn_id in rxn_ids:
            try:
                rxn = model.reactions.get_by_id(rxn_id)
                genes = [gene.id for gene in rxn.genes]
                reaction_equation = rxn.reaction
            except Exception:
                genes = []
                reaction_equation = ""
            row = [rxn_id, reaction_equation, genes] + [df_new.loc[rxn_id, col] for col in flux_cols]
            data.append(row)
        return data
    
    upregulated_data = build_output_data(upregulated)
    downregulated_data = build_output_data(downregulated)
    knockout_data = build_output_data(knockout)
    
    # 创建DataFrame并输出到Excel
    columns = ['Reaction ID', 'Reaction', 'Genes'] + [f'Flux {i}' for i in flux_cols]
    upregulated_df = pd.DataFrame(upregulated_data, columns=columns)
    downregulated_df = pd.DataFrame(downregulated_data, columns=columns)
    knockout_df = pd.DataFrame(knockout_data, columns=columns)
    
    output_filename = f'FSEOF_{objective_reaction}_results.xlsx'
    with pd.ExcelWriter(output_filename) as writer:
        upregulated_df.to_excel(writer, sheet_name='Upregulated', index=False)
        downregulated_df.to_excel(writer, sheet_name='Downregulated', index=False)
        knockout_df.to_excel(writer, sheet_name='Knockout', index=False)
    
    print(f"FSEOF分析完成，结果已保存到 {output_filename}")

if __name__ == "__main__":
    # 命令行参数解析
    parser = argparse.ArgumentParser(description='进行FSEOF分析')
    parser.add_argument('model_path', type=str, help='SBML模型文件路径')
    parser.add_argument('objective_reaction', type=str, help='目标反应ID')
    parser.add_argument('biomass_reaction', type=str, help='生物量反应ID')
    args = parser.parse_args()
    
    # 调用FSEOF分析函数
    fseof_analysis(args.model_path, args.objective_reaction, args.biomass_reaction)
