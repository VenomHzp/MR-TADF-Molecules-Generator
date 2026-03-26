import csv
import os
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor

# 核心步骤：强制启用薛定谔捐赠的高级坐标生成引擎 (CoordGen)
# 这对避免大分子、多环和长侧链的重叠有奇效
rdDepictor.SetPreferCoordGen(True)

def visualize_smiles_from_csv(csv_file, output_dir="output_images", chunk_size=50):
    """
    读取 CSV 文件中的 SMILES，每 50 个生成一张高清大图。
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    smiles_list = []
    ids_list = []
    
    # 1. 读取 CSV
    print(f"正在读取 {csv_file}...")
    with open(csv_file, mode='r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            smiles_list.append(row["SMILES"])
            ids_list.append(row["ID"])
            
    total_mols = len(smiles_list)
    print(f"共读取到 {total_mols} 个分子。开始渲染...")

    # 2. 分块处理 (每 chunk_size 个分子为一组)
    for i in range(0, total_mols, chunk_size):
        chunk_smiles = smiles_list[i : i + chunk_size]
        chunk_ids = ids_list[i : i + chunk_size]
        
        mols = []
        valid_ids = []
        for sm, mol_id in zip(chunk_smiles, chunk_ids):
            mol = Chem.MolFromSmiles(sm)
            if mol:
                mols.append(mol)
                valid_ids.append(mol_id)
        
        # 3. 绘图设置
        # molsPerRow=5: 每行 5 个分子
        # subImgSize=(500, 500): 大幅度拉大每个分子的独立画布尺寸 (默认通常是 200x200)
        batch_num = (i // chunk_size) + 1
        output_path = os.path.join(output_dir, f"mols_batch_{batch_num}.png")
        
        try:
            img = Draw.MolsToGridImage(
                mols,
                molsPerRow=5,
                subImgSize=(500, 500), 
                legends=valid_ids,
                useSVG=False # 设为 True 会保存为无损矢量图，设为 False 为 PNG
            )
            # 保存图像
            img.save(output_path)
            print(f"已保存图像: {output_path} (包含 {len(mols)} 个分子)")
        except Exception as e:
            print(f"保存第 {batch_num} 批图像时出错: {e}")

if __name__ == "__main__":
    csv_filename = "novel_mr_skeletons.csv"
    if os.path.exists(csv_filename):
        visualize_smiles_from_csv(csv_filename, chunk_size=50)
    else:
        print(f"找不到文件 {csv_filename}，请先运行主生成脚本。")