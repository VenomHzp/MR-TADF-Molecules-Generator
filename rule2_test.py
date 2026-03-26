from rdkit import Chem
import random

class AdvancedMRMutator:
    def __init__(self):
        # 法则 2 的候选桥连基团及其权重
        self.bridge_weights = {
            'N-Ph': 0.30,
            'B-Ph': 0.30,
            'O': 0.066,
            'S': 0.066,
            'Se': 0.066,
            'O=S=O (砜基)': 0.066,
            'P=O (磷氧)': 0.066,
            'C=O (羰基)': 0.066
        }

    # --- 辅助方法：安全添加苯环 (采用 Kekule 模式防止芳香性崩溃) ---
    def _add_phenyl_ring(self, rwmol, attach_idx):
        start_idx = rwmol.GetNumAtoms()
        # 添加 6 个碳原子
        for i in range(6):
            rwmol.AddAtom(Chem.Atom(6))
        # 连接到主骨架
        rwmol.AddBond(attach_idx, start_idx, Chem.BondType.SINGLE)
        # 用交替的单双键（Kekule式）闭合苯环，后续 Sanitize 会自动识别为芳香环
        rwmol.AddBond(start_idx, start_idx+1, Chem.BondType.DOUBLE)
        rwmol.AddBond(start_idx+1, start_idx+2, Chem.BondType.SINGLE)
        rwmol.AddBond(start_idx+2, start_idx+3, Chem.BondType.DOUBLE)
        rwmol.AddBond(start_idx+3, start_idx+4, Chem.BondType.SINGLE)
        rwmol.AddBond(start_idx+4, start_idx+5, Chem.BondType.DOUBLE)
        rwmol.AddBond(start_idx+5, start_idx, Chem.BondType.SINGLE)

    # --- 法则 2：邻接桥连平面化 (升级版：纯拓扑寻路) ---
    def rule2_planarization(self, mol):
        """
        寻找与中心杂原子间隔 1 个原子的苯环 C-H 位点。
        使用“纯拓扑寻路 + 同环排斥锁”精确定位海湾区，彻底抛弃不稳定的 3D 测距。
        """
        candidates = []
        target_heteros = [5, 7, 8, 16, 34, 15] # B, N, O, S, Se, P

        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() in target_heteros:
                # 寻找间隔 1 个原子的苯环 ipso-碳
                ipsos = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIsAromatic()]
                
                if len(ipsos) >= 2:
                    for i in range(len(ipsos)):
                        for j in range(i+1, len(ipsos)):
                            ipso1, ipso2 = ipsos[i], ipsos[j]

                            # ==================================================
                            # 终极护城河 1：同环排斥锁 (Fused Ring Lock)
                            # 如果这三个原子已经在一个刚性的 5、6、7 元环里，
                            # 说明它们的外侧键必定是发散的（V 型），直接跳过！
                            # ==================================================
                            in_same_ring = False
                            for ring in mol.GetRingInfo().AtomRings():
                                if atom.GetIdx() in ring and ipso1.GetIdx() in ring and ipso2.GetIdx() in ring:
                                    if len(ring) <= 7:
                                        in_same_ring = True
                                        break
                            if in_same_ring:
                                continue

                            # 寻找对应的 ortho-碳（必须有 H 位点）
                            orthos1 = [n for n in ipso1.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIsAromatic() and n.GetTotalNumHs() > 0 and n.GetIdx() != atom.GetIdx()]
                            orthos2 = [n for n in ipso2.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIsAromatic() and n.GetTotalNumHs() > 0 and n.GetIdx() != atom.GetIdx()]

                            # 终极护城河 2：对称去重锁
                            bay_found_for_this_pair = False
                            for o1 in orthos1:
                                for o2 in orthos2:
                                    if bay_found_for_this_pair:
                                        break
                                        
                                    # 纯拓扑寻路：最短路径必须精确等于 5 个原子，且必须穿过中心杂原子
                                    path = Chem.rdmolops.GetShortestPath(mol, o1.GetIdx(), o2.GetIdx())
                                    if len(path) == 5 and atom.GetIdx() in path:
                                        candidates.append((o1.GetIdx(), o2.GetIdx()))
                                        bay_found_for_this_pair = True

        if not candidates:
            return None, "未找到符合条件的纯拓扑海湾区 C-H 位点"

        # 随机选择一对候选 C 原子
        c1_idx, c2_idx = random.choice(candidates)
        rwmol = Chem.RWMol(mol)

        # 按概率抽取桥连片段
        bridge_types = list(self.bridge_weights.keys())
        weights = list(self.bridge_weights.values())
        selected_bridge = random.choices(bridge_types, weights=weights, k=1)[0]

        # ========= 执行复杂的片段桥连手术 =========
        if selected_bridge == 'O':
            idx = rwmol.AddAtom(Chem.Atom(8))
            rwmol.AddBond(c1_idx, idx, Chem.BondType.SINGLE)
            rwmol.AddBond(c2_idx, idx, Chem.BondType.SINGLE)
            
        elif selected_bridge == 'S':
            idx = rwmol.AddAtom(Chem.Atom(16))
            rwmol.AddBond(c1_idx, idx, Chem.BondType.SINGLE)
            rwmol.AddBond(c2_idx, idx, Chem.BondType.SINGLE)
            
        elif selected_bridge == 'Se':
            idx = rwmol.AddAtom(Chem.Atom(34))
            rwmol.AddBond(c1_idx, idx, Chem.BondType.SINGLE)
            rwmol.AddBond(c2_idx, idx, Chem.BondType.SINGLE)
            
        elif selected_bridge == 'C=O (羰基)':
            c_idx = rwmol.AddAtom(Chem.Atom(6))
            o_idx = rwmol.AddAtom(Chem.Atom(8))
            rwmol.AddBond(c1_idx, c_idx, Chem.BondType.SINGLE)
            rwmol.AddBond(c2_idx, c_idx, Chem.BondType.SINGLE)
            rwmol.AddBond(c_idx, o_idx, Chem.BondType.DOUBLE)
            
        elif selected_bridge == 'O=S=O (砜基)':
            s_idx = rwmol.AddAtom(Chem.Atom(16))
            o1 = rwmol.AddAtom(Chem.Atom(8))
            o2 = rwmol.AddAtom(Chem.Atom(8))
            rwmol.AddBond(c1_idx, s_idx, Chem.BondType.SINGLE)
            rwmol.AddBond(c2_idx, s_idx, Chem.BondType.SINGLE)
            rwmol.AddBond(s_idx, o1, Chem.BondType.DOUBLE)
            rwmol.AddBond(s_idx, o2, Chem.BondType.DOUBLE)
            
        elif selected_bridge == 'P=O (磷氧)':
            p_idx = rwmol.AddAtom(Chem.Atom(15))
            o_idx = rwmol.AddAtom(Chem.Atom(8))
            rwmol.AddBond(c1_idx, p_idx, Chem.BondType.SINGLE)
            rwmol.AddBond(c2_idx, p_idx, Chem.BondType.SINGLE)
            rwmol.AddBond(p_idx, o_idx, Chem.BondType.DOUBLE)
            self._add_phenyl_ring(rwmol, p_idx) # 给 P 加上苯环满足价态
            
        elif selected_bridge == 'N-Ph':
            n_idx = rwmol.AddAtom(Chem.Atom(7))
            rwmol.AddBond(c1_idx, n_idx, Chem.BondType.SINGLE)
            rwmol.AddBond(c2_idx, n_idx, Chem.BondType.SINGLE)
            self._add_phenyl_ring(rwmol, n_idx)
            
        elif selected_bridge == 'B-Ph':
            b_idx = rwmol.AddAtom(Chem.Atom(5))
            rwmol.AddBond(c1_idx, b_idx, Chem.BondType.SINGLE)
            rwmol.AddBond(c2_idx, b_idx, Chem.BondType.SINGLE)
            self._add_phenyl_ring(rwmol, b_idx)

        # 校验化合价和芳香性
        try:
            Chem.SanitizeMol(rwmol)
            return rwmol.GetMol(), selected_bridge
        except Exception:
            return None, selected_bridge

# ==========================================
# 运行测试：DABNA 骨架的完全锁死
# ==========================================
if __name__ == "__main__":
    mutator = AdvancedMRMutator()
    
    dabna_smiles = "C1=CC(=CC=C1)N1C2C(B3C4=CC=CC=C4N(C4C=CC=CC=4)C4C=CC=C1C=43)=CC=CC=2"
    mol = Chem.MolFromSmiles(dabna_smiles)
    
    print(f"【初始状态】加载经典 DABNA 骨架: {dabna_smiles}")
    print("-" * 50)
    
    # 连续运行 5 次，展示插入不同基团的效果
    for i in range(5):
        new_mol, bridge_used = mutator.rule2_planarization(mol)
        if new_mol:
            new_smiles = Chem.MolToSmiles(new_mol)
            print(f"[测试 {i+1}] 成功捕获海湾区！")
            print(f"   插入片段: {bridge_used}")
            print(f"   生成新型全刚性 MR 骨架: {new_smiles}\n")
        else:
            print(f"[测试 {i+1}] 变异失败，尝试插入 {bridge_used} 时化合价未通过\n")