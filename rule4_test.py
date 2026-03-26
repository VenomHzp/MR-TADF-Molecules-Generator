from rdkit import Chem

def apply_rule4_to_all_bays(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return "无效的 SMILES"

    target_heteros = [5, 7, 8, 16, 34, 15] 
    bonds_to_form = set()
    
    # 【核心升级】：建立全局消耗池。一个碳原子只能结一次婚！
    global_used_carbons = set() 
    
    print(f"\n🚀 开始执行全海湾闭合手术 (法则 4 的极限测试)...")
    
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in target_heteros:
            ipsos = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIsAromatic()]
            if len(ipsos) >= 2:
                for i in range(len(ipsos)):
                    for j in range(i+1, len(ipsos)):
                        ipso1, ipso2 = ipsos[i], ipsos[j]
                        
                        # 同环排斥锁 (防 V 型)
                        in_same_ring = False
                        for ring in mol.GetRingInfo().AtomRings():
                            if atom.GetIdx() in ring and ipso1.GetIdx() in ring and ipso2.GetIdx() in ring:
                                if len(ring) <= 7:
                                    in_same_ring = True
                                    break
                        if in_same_ring:
                            continue

                        orthos1 = [n for n in ipso1.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIsAromatic() and n.GetTotalNumHs() > 0 and n.GetIdx() != atom.GetIdx()]
                        orthos2 = [n for n in ipso2.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIsAromatic() and n.GetTotalNumHs() > 0 and n.GetIdx() != atom.GetIdx()]

                        bay_found_for_this_pair = False
                        for o1 in orthos1:
                            for o2 in orthos2:
                                if bay_found_for_this_pair:
                                    break
                                
                                # ==================================================
                                # 终极护城河：消耗锁 (Consumption Lock)
                                # 彻底杜绝同一个碳原子被劈腿接两次！
                                # ==================================================
                                if o1.GetIdx() in global_used_carbons or o2.GetIdx() in global_used_carbons:
                                    continue # 已经被别的海湾用掉了，跳过！
                                    
                                path = Chem.rdmolops.GetShortestPath(mol, o1.GetIdx(), o2.GetIdx())
                                if len(path) == 5 and atom.GetIdx() in path:
                                    # 记录要连的键
                                    bond_pair = tuple(sorted((o1.GetIdx(), o2.GetIdx())))
                                    bonds_to_form.add(bond_pair)
                                    
                                    # 将这两个碳扔进消耗池，宣布它们名花有主
                                    global_used_carbons.add(o1.GetIdx())
                                    global_used_carbons.add(o2.GetIdx())
                                    
                                    bay_found_for_this_pair = True

    if not bonds_to_form:
        return "未发现任何海湾，分子已是全刚性。"

    print(f"🎯 扫描完毕！共锁定 {len(bonds_to_form)} 个真实海湾，准备闭环：")
    
    rwmol = Chem.RWMol(mol)
    for c1, c2 in bonds_to_form:
        print(f"   -> 正在焊接碳原子 [{c1}] 和 [{c2}] ...")
        rwmol.AddBond(c1, c2, Chem.BondType.SINGLE)
        
    print("🔄 正在交由 RDKit 进行全分子芳香性重排与审查...")
    try:
        Chem.SanitizeMol(rwmol)
        final_smiles = Chem.MolToSmiles(rwmol.GetMol())
        print("✅ 审查通过！生成了完美的超级刚性稠环体系。")
        return final_smiles
    except Exception as e:
        return f"❌ 闭环失败！物理张力过大或共轭电子数无法满足芳香性: {e}"

if __name__ == "__main__":
    dabna_smiles = "C1=CC(=CC=C1)N1C2C(B3C4=CC=CC=C4N(C4C=CC=CC=4)C4C=CC=C1C=43)=CC=CC=2"
    
    final_smi = apply_rule4_to_all_bays(dabna_smiles)
    
    print("\n" + "="*60)
    print("【初始骨架 SMILES】:", dabna_smiles)
    print("【极限闭环 SMILES】:", final_smi)
    print("="*60)