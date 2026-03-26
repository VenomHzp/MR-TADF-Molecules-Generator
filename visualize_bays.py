from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
import os

rdDepictor.SetPreferCoordGen(True)

def visualize_bay_regions_ultimate(smiles, output_filename="bay_check_ultimate.png"):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        print("无效的 SMILES")
        return

    rdDepictor.Compute2DCoords(mol)

    target_heteros = [5, 7, 8, 16, 34, 15] 
    bay_carbons = set()
    anchor_heteros = set()
    found_pairs = []

    print(f"\n🔍 开始扫描分子 (启用纯图论 + 同环排斥锁)...")
    
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in target_heteros:
            ipsos = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIsAromatic()]
            if len(ipsos) >= 2:
                for i in range(len(ipsos)):
                    for j in range(i+1, len(ipsos)):
                        ipso1, ipso2 = ipsos[i], ipsos[j]
                        
                        # ==================================================
                        # 终极护城河：同环排斥锁 (Fused Ring Lock)
                        # 如果这三个原子已经在一个刚性的 5、6、7 元环里，
                        # 说明它们的外侧键必定是发散的（V 型），直接跳过这对苯环！
                        # ==================================================
                        in_same_ring = False
                        for ring in mol.GetRingInfo().AtomRings():
                            if atom.GetIdx() in ring and ipso1.GetIdx() in ring and ipso2.GetIdx() in ring:
                                if len(ring) <= 7:
                                    in_same_ring = True
                                    break
                        if in_same_ring:
                            continue # 直接拦截假阳性！

                        # 继续寻找这对未被底层骨架锁死的苯环的海湾
                        orthos1 = [n for n in ipso1.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIsAromatic() and n.GetTotalNumHs() > 0 and n.GetIdx() != atom.GetIdx()]
                        orthos2 = [n for n in ipso2.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIsAromatic() and n.GetTotalNumHs() > 0 and n.GetIdx() != atom.GetIdx()]

                        bay_found_for_this_pair = False

                        for o1 in orthos1:
                            for o2 in orthos2:
                                if bay_found_for_this_pair:
                                    break
                                    
                                path = Chem.rdmolops.GetShortestPath(mol, o1.GetIdx(), o2.GetIdx())
                                
                                if len(path) == 5 and atom.GetIdx() in path:
                                    bay_carbons.add(o1.GetIdx())
                                    bay_carbons.add(o2.GetIdx())
                                    anchor_heteros.add(atom.GetIdx())
                                    found_pairs.append((atom.GetSymbol(), o1.GetIdx(), o2.GetIdx()))
                                    
                                    bay_found_for_this_pair = True 

    if not bay_carbons:
        print("✅ 未找到海湾区")
        return

    print(f"🎯 扫描完毕！算法精准锁定 {len(found_pairs)} 对真实海湾碳原子：")
    for pair in found_pairs:
        print(f"   - 中心 [{pair[0]}] 锚点 -> 碳原子 {pair[1]} 与 {pair[2]}")

    # ==========================================
    # 渲染带原子序号的高清图
    # ==========================================
    highlight_atoms = list(bay_carbons) + list(anchor_heteros)
    highlight_colors = {idx: (1.0, 0.4, 0.4) for idx in bay_carbons} 
    highlight_colors.update({idx: (0.4, 0.6, 1.0) for idx in anchor_heteros}) 

    d2d = Draw.rdMolDraw2D.MolDraw2DCairo(600, 600)
    opts = d2d.drawOptions()
    opts.addAtomIndices = True     
    opts.clearBackground = True
    opts.highlightBondWidthMultiplier = 12 

    Draw.rdMolDraw2D.PrepareAndDrawMolecule(
        d2d, mol, highlightAtoms=highlight_atoms, highlightAtomColors=highlight_colors
    )
    d2d.FinishDrawing()
    
    with open(output_filename, 'wb') as f:
        f.write(d2d.GetDrawingText())
        
    print(f"🖼️ 可视化结果已保存至: {os.path.abspath(output_filename)}\n")

if __name__ == "__main__":
    dabna_smiles = "C1=CC(=CC=C1)N1C2C(B3C4=CC=CC=C4N(C4C=CC=CC=4)C4C=CC=C1C=43)=CC=CC=2"
    visualize_bay_regions_ultimate(dabna_smiles, "dabna_bay_ultimate.png")