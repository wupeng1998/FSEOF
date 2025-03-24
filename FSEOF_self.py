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
    initial_flux = 0
    number_of_results = 10
    levels = [initial_flux + (i + 1) * (max_flux - initial_flux) / number_of_results for i in range(number_of_results)]
    print(f"通量水平: {levels}")
    
    # 计算不同通量水平下的代谢网络
    df_list = []
    for f in levels:
        with model:
            model.reactions.get_by_id(objective_reaction).bounds = (f, f)
            model.objective = biomass_reaction_obj
            pfba_solution = pfba(model)
            pfba_solution = pfba_solution.to_frame()
            pfba_solution.drop('reduced_costs', axis=1, inplace=True)
            df_list.append(pfba_solution)
    
    df_new = pd.concat(df_list, axis=1)
    df_new.columns = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    df_new.fillna(value=0, inplace=True)
    
    # 数据处理
    remove_rea = ['GLCtex_copy1', 'NH4tex', 'NH4tpp', 'H2Otex', 'H2Otpp', 'EX_h2o_e', 'O2tpp', 'O2tex', 'EX_o2_e', 'CO2tpp', 'CO2tex', 'Htex', 'EX_h_e', 'ADD_H_c-tex', 'ADD_EX_h_c']
    weight = ['ATPM', 'TPI', 'ENO', 'PGM', 'PGK', 'GLCtex_copy1', 'NH4tex', 'NH4tpp', 'EX_nh4_e', 'H2Otex', 'H2Otpp', 'EX_h2o_e', 'O2tpp', 'O2tex', 'EX_o2_e', 'CO2tpp', 'CO2tex', 'EX_co2_e', 'Htex', 'EX_h_e', 'ADD_H_c-tex', 'ADD_EX_h_c']
    lose_weight = ['F6PA', 'FBA3', 'EDA', 'XYLI2']
    
    # 去掉通量绝对值小于等于1e-5的行
    for i in df_new.index:
        n = 0
        for v in df_new.loc[i]:
            if abs(v) <= 1e-5:
                n += 1
        if n == 10:
            df_new.drop(i, inplace=True)
    
        for re_rea in remove_rea:      #删除指定的反应行
            if i == re_rea:
                try:
                    df_new.drop(i,inplace=True)
                except:
                    pass
    
    # 调整权重
    for ii in df_new.index:
        for weight_1 in weight:
            if ii == weight_1:
                df_new.loc[ii, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]] = df_new.loc[ii, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]] * 1000
        for lose_weight_1 in lose_weight:
            if ii == lose_weight_1:
                df_new.loc[ii, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]] = df_new.loc[ii, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]] / 50
    
    df_new = df_new.round(5)
    
    # 插入反应方程式列
    fseof_rea = []
    for i in df_new.index:
        r = model.reactions.get_by_id(i)
        fseof_rea.append(r.reaction)
    df_new.insert(0, 'reactions', fseof_rea)
    
    # 筛选强化靶点
    upregulated = []
    for key_up, row_up in df_new.iterrows():
        max_value = max(list(abs(row_up[1:11])))
        value_initial = df_new.loc[key_up][1]
        value_end = df_new.loc[key_up][10]
        
        if value_initial * value_end >= 0 and abs(value_end) > abs(value_initial):
            if max_value == abs(value_end):
                upregulated.append(key_up)
            elif ((max_value - abs(value_end)) / max_value) < 0.1:
                upregulated.append(key_up)
    
    # 筛选弱化靶点
    downregulated = []
    for key_down, row_down in df_new.iterrows():
        min_value = min(list(abs(row_down[1:11])))
        value_initial = df_new.loc[key_down][1]
        value_end = df_new.loc[key_down][10]
        
        if value_initial * value_end >= 0 and abs(value_end) < abs(value_initial):
            if min_value == abs(value_end):
                downregulated.append(key_down)
    
    # 筛选敲除靶点
    knockout = []
    for key_ko, row_ko in df_new.iterrows():
        value_end = df_new.loc[key_ko][10]
        if value_end == 0 and df_new.loc[key_ko][1] != 0:
            knockout.append(key_ko)
    
    # 准备输出数据
    upregulated_data = []
    for i in upregulated:
        reaction = model.reactions.get_by_id(i)
        row = [
            i,
            reaction.reaction,
            [ge.id for ge in reaction.genes],
            df_new.loc[i][1],
            df_new.loc[i][2],
            df_new.loc[i][3],
            df_new.loc[i][4],
            df_new.loc[i][5],
            df_new.loc[i][6],
            df_new.loc[i][7],
            df_new.loc[i][8],
            df_new.loc[i][9],
            df_new.loc[i][10]
        ]
        upregulated_data.append(row)
    
    downregulated_data = []
    for i in downregulated:
        reaction = model.reactions.get_by_id(i)
        row = [
            i,
            reaction.reaction,
            [ge.id for ge in reaction.genes],
            df_new.loc[i][1],
            df_new.loc[i][2],
            df_new.loc[i][3],
            df_new.loc[i][4],
            df_new.loc[i][5],
            df_new.loc[i][6],
            df_new.loc[i][7],
            df_new.loc[i][8],
            df_new.loc[i][9],
            df_new.loc[i][10]
        ]
        downregulated_data.append(row)
    
    knockout_data = []
    for i in knockout:
        reaction = model.reactions.get_by_id(i)
        row = [
            i,
            reaction.reaction,
            [ge.id for ge in reaction.genes],
            df_new.loc[i][1],
            df_new.loc[i][2],
            df_new.loc[i][3],
            df_new.loc[i][4],
            df_new.loc[i][5],
            df_new.loc[i][6],
            df_new.loc[i][7],
            df_new.loc[i][8],
            df_new.loc[i][9],
            df_new.loc[i][10]
        ]
        knockout_data.append(row)
    
    # 创建DataFrame
    upregulated_df = pd.DataFrame(upregulated_data, columns=[
        'Reaction ID', 'Reaction', 'Genes', 'Flux 1', 'Flux 2', 'Flux 3', 'Flux 4', 'Flux 5',
        'Flux 6', 'Flux 7', 'Flux 8', 'Flux 9', 'Flux 10'
    ])
    
    downregulated_df = pd.DataFrame(downregulated_data, columns=[
        'Reaction ID', 'Reaction', 'Genes', 'Flux 1', 'Flux 2', 'Flux 3', 'Flux 4', 'Flux 5',
        'Flux 6', 'Flux 7', 'Flux 8', 'Flux 9', 'Flux 10'
    ])
    
    knockout_df = pd.DataFrame(knockout_data, columns=[
        'Reaction ID', 'Reaction', 'Genes', 'Flux 1', 'Flux 2', 'Flux 3', 'Flux 4', 'Flux 5',
        'Flux 6', 'Flux 7', 'Flux 8', 'Flux 9', 'Flux 10'
    ])
    
    # 输出到Excel
    with pd.ExcelWriter(f'FSEOF_{objective_reaction}_results.xlsx') as writer:
        upregulated_df.to_excel(writer, sheet_name='Upregulated', index=False)
        downregulated_df.to_excel(writer, sheet_name='Downregulated', index=False)
        knockout_df.to_excel(writer, sheet_name='Knockout', index=False)
    
    print(f"FSEOF分析完成，结果已保存到 FSEOF_{objective_reaction}_results.xlsx")

if __name__ == "__main__":
    # 命令行参数解析
    parser = argparse.ArgumentParser(description='进行FSEOF分析')
    parser.add_argument('model_path', type=str, help='SBML模型文件路径')
    parser.add_argument('objective_reaction', type=str, help='目标反应ID')
    parser.add_argument('biomass_reaction', type=str, help='生物量反应ID')
    args = parser.parse_args()
    
    # 调用FSEOF分析函数
    fseof_analysis(args.model_path, args.objective_reaction, args.biomass_reaction)