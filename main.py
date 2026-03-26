import os
import csv
import time
import random
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor

# 强制开启薛定谔高级 2D 坐标生成引擎
rdDepictor.SetPreferCoordGen(True)
# 屏蔽 RDKit 底层刷屏警告，保持终端清爽
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

class MRMutationEngine:
    def __init__(self):
        # 法则 2 的候选桥连基团
        self.bridge_weights = {
            'N-Ph': 0.30, 'B-Ph': 0.30, 'O': 0.066, 'S': 0.066,
            'Se': 0.066, 'O=S=O': 0.066, 'P=O': 0.066, 'C=O': 0.066
        }
        # 法则 3 的候选扩增片段
        self.frag_types = ['Benzene', 'Furan', 'Benzofuran', 'N-Ph-Indole']
        self.frag_weights = [0.55, 0.15, 0.15, 0.15]

    # ==========================================
    # 辅助工具包：无损拼接与缝合
    # ==========================================
    def _add_phenyl_ring(self, rwmol, attach_idx):
        start_idx = rwmol.GetNumAtoms()
        for i in range(6): rwmol.AddAtom(Chem.Atom(6))
        rwmol.AddBond(attach_idx, start_idx, Chem.BondType.SINGLE)
        rwmol.AddBond(start_idx, start_idx+1, Chem.BondType.DOUBLE)
        rwmol.AddBond(start_idx+1, start_idx+2, Chem.BondType.SINGLE)
        rwmol.AddBond(start_idx+2, start_idx+3, Chem.BondType.DOUBLE)
        rwmol.AddBond(start_idx+3, start_idx+4, Chem.BondType.SINGLE)
        rwmol.AddBond(start_idx+4, start_idx+5, Chem.BondType.DOUBLE)
        rwmol.AddBond(start_idx+5, start_idx, Chem.BondType.SINGLE)

    def _fuse_benzene(self, rwmol, a1, a2):
        new_atoms = [rwmol.AddAtom(Chem.Atom(6)) for _ in range(4)]
        for idx in new_atoms: rwmol.GetAtomWithIdx(idx).SetIsAromatic(True)
        rwmol.AddBond(a1, new_atoms[0], Chem.BondType.AROMATIC)
        rwmol.AddBond(new_atoms[0], new_atoms[1], Chem.BondType.AROMATIC)
        rwmol.AddBond(new_atoms[1], new_atoms[2], Chem.BondType.AROMATIC)
        rwmol.AddBond(new_atoms[2], new_atoms[3], Chem.BondType.AROMATIC)
        rwmol.AddBond(new_atoms[3], a2, Chem.BondType.AROMATIC)

    def _fuse_furan(self, rwmol, a1, a2):
        o_idx, c1_idx, c2_idx = rwmol.AddAtom(Chem.Atom(8)), rwmol.AddAtom(Chem.Atom(6)), rwmol.AddAtom(Chem.Atom(6))
        for idx in [o_idx, c1_idx, c2_idx]: rwmol.GetAtomWithIdx(idx).SetIsAromatic(True)
        if random.choice([True, False]):
            rwmol.AddBond(a1, o_idx, Chem.BondType.AROMATIC)
            rwmol.AddBond(o_idx, c1_idx, Chem.BondType.AROMATIC)
            rwmol.AddBond(c1_idx, c2_idx, Chem.BondType.AROMATIC)
            rwmol.AddBond(c2_idx, a2, Chem.BondType.AROMATIC)
        else:
            rwmol.AddBond(a2, o_idx, Chem.BondType.AROMATIC)
            rwmol.AddBond(o_idx, c1_idx, Chem.BondType.AROMATIC)
            rwmol.AddBond(c1_idx, c2_idx, Chem.BondType.AROMATIC)
            rwmol.AddBond(c2_idx, a1, Chem.BondType.AROMATIC)

    def _fuse_benzofuran(self, rwmol, a1, a2):
        o_idx = rwmol.AddAtom(Chem.Atom(8))
        rwmol.GetAtomWithIdx(o_idx).SetIsAromatic(True)
        bc = [rwmol.AddAtom(Chem.Atom(6)) for _ in range(6)]
        for idx in bc: rwmol.GetAtomWithIdx(idx).SetIsAromatic(True)
        for i in range(6): rwmol.AddBond(bc[i], bc[(i+1)%6], Chem.BondType.AROMATIC)
        if random.choice([True, False]):
            rwmol.AddBond(a1, o_idx, Chem.BondType.AROMATIC)
            rwmol.AddBond(o_idx, bc[0], Chem.BondType.AROMATIC)
            rwmol.AddBond(bc[1], a2, Chem.BondType.AROMATIC)
        else:
            rwmol.AddBond(a2, o_idx, Chem.BondType.AROMATIC)
            rwmol.AddBond(o_idx, bc[0], Chem.BondType.AROMATIC)
            rwmol.AddBond(bc[1], a1, Chem.BondType.AROMATIC)

    def _fuse_nphenylindole(self, rwmol, a1, a2):
        n_idx = rwmol.AddAtom(Chem.Atom(7))
        rwmol.GetAtomWithIdx(n_idx).SetIsAromatic(True)
        bc, ph = [rwmol.AddAtom(Chem.Atom(6)) for _ in range(6)], [rwmol.AddAtom(Chem.Atom(6)) for _ in range(6)]
        for idx in bc + ph: rwmol.GetAtomWithIdx(idx).SetIsAromatic(True)
        for i in range(6):
            rwmol.AddBond(bc[i], bc[(i+1)%6], Chem.BondType.AROMATIC)
            rwmol.AddBond(ph[i], ph[(i+1)%6], Chem.BondType.AROMATIC)
        rwmol.AddBond(n_idx, ph[0], Chem.BondType.SINGLE)
        if random.choice([True, False]):
            rwmol.AddBond(a1, n_idx, Chem.BondType.AROMATIC)
            rwmol.AddBond(n_idx, bc[0], Chem.BondType.AROMATIC)
            rwmol.AddBond(bc[1], a2, Chem.BondType.AROMATIC)
        else:
            rwmol.AddBond(a2, n_idx, Chem.BondType.AROMATIC)
            rwmol.AddBond(n_idx, bc[0], Chem.BondType.AROMATIC)
            rwmol.AddBond(bc[1], a1, Chem.BondType.AROMATIC)

    # ==========================================
    # 四大变异法则 (带全局安全锁)
    # ==========================================
    def rule1_n_doping(self, mol):
        rwmol = Chem.RWMol(mol)
        c_cands = sorted([a.GetIdx() for a in rwmol.GetAtoms() if a.GetAtomicNum() == 6 and a.GetIsAromatic()])
        safe_cands = []
        for c in c_cands:
            atom = rwmol.GetAtomWithIdx(c)
            if atom.GetTotalNumHs() == 0: continue # 保护桥头碳
            neighbors = [n.GetAtomicNum() for n in atom.GetNeighbors()]
            if all(num in [6, 1] for num in neighbors): # 绝对杂原子隔离
                safe_cands.append(c)
        if not safe_cands: return mol, False, "无合法碳"
        
        target = random.choice(safe_cands)
        rwmol.GetAtomWithIdx(target).SetAtomicNum(7)
        try:
            Chem.SanitizeMol(rwmol)
            return rwmol.GetMol(), True, "Rule1"
        except: return mol, False, "化合价错误"

    def rule2_planarization(self, mol):
        candidates, target_heteros = [], [5, 7, 8, 16, 34, 15]
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() in target_heteros:
                ipsos = sorted([n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIsAromatic()], key=lambda x: x.GetIdx())
                if len(ipsos) >= 2:
                    for i in range(len(ipsos)):
                        for j in range(i+1, len(ipsos)):
                            ip1, ip2 = ipsos[i], ipsos[j]
                            
                            # 同环排斥锁
                            in_same_ring = False
                            for ring in mol.GetRingInfo().AtomRings():
                                if atom.GetIdx() in ring and ip1.GetIdx() in ring and ip2.GetIdx() in ring:
                                    if len(ring) <= 7: in_same_ring = True; break
                            if in_same_ring: continue

                            orthos1 = [n for n in ip1.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIsAromatic() and n.GetTotalNumHs() > 0 and n.GetIdx() != atom.GetIdx()]
                            orthos2 = [n for n in ip2.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIsAromatic() and n.GetTotalNumHs() > 0 and n.GetIdx() != atom.GetIdx()]

                            bay_found = False
                            for o1 in orthos1:
                                for o2 in orthos2:
                                    if bay_found: break
                                    path = Chem.rdmolops.GetShortestPath(mol, o1.GetIdx(), o2.GetIdx())
                                    if len(path) == 5 and atom.GetIdx() in path:
                                        candidates.append((o1.GetIdx(), o2.GetIdx()))
                                        bay_found = True
        if not candidates: return mol, False, "无海湾"

        c1_idx, c2_idx = random.choice(candidates)
        rwmol = Chem.RWMol(mol)
        bridge = random.choices(list(self.bridge_weights.keys()), weights=list(self.bridge_weights.values()), k=1)[0]
        
        if bridge == 'O':
            idx = rwmol.AddAtom(Chem.Atom(8))
            rwmol.AddBond(c1_idx, idx, Chem.BondType.SINGLE); rwmol.AddBond(c2_idx, idx, Chem.BondType.SINGLE)
        elif bridge == 'S':
            idx = rwmol.AddAtom(Chem.Atom(16))
            rwmol.AddBond(c1_idx, idx, Chem.BondType.SINGLE); rwmol.AddBond(c2_idx, idx, Chem.BondType.SINGLE)
        elif bridge == 'Se':
            idx = rwmol.AddAtom(Chem.Atom(34))
            rwmol.AddBond(c1_idx, idx, Chem.BondType.SINGLE); rwmol.AddBond(c2_idx, idx, Chem.BondType.SINGLE)
        elif bridge == 'C=O':
            c_idx, o_idx = rwmol.AddAtom(Chem.Atom(6)), rwmol.AddAtom(Chem.Atom(8))
            rwmol.AddBond(c1_idx, c_idx, Chem.BondType.SINGLE); rwmol.AddBond(c2_idx, c_idx, Chem.BondType.SINGLE)
            rwmol.AddBond(c_idx, o_idx, Chem.BondType.DOUBLE)
        elif bridge == 'O=S=O':
            s_idx, o1, o2 = rwmol.AddAtom(Chem.Atom(16)), rwmol.AddAtom(Chem.Atom(8)), rwmol.AddAtom(Chem.Atom(8))
            rwmol.AddBond(c1_idx, s_idx, Chem.BondType.SINGLE); rwmol.AddBond(c2_idx, s_idx, Chem.BondType.SINGLE)
            rwmol.AddBond(s_idx, o1, Chem.BondType.DOUBLE); rwmol.AddBond(s_idx, o2, Chem.BondType.DOUBLE)
        elif bridge == 'P=O':
            p_idx, o_idx = rwmol.AddAtom(Chem.Atom(15)), rwmol.AddAtom(Chem.Atom(8))
            rwmol.AddBond(c1_idx, p_idx, Chem.BondType.SINGLE); rwmol.AddBond(c2_idx, p_idx, Chem.BondType.SINGLE)
            rwmol.AddBond(p_idx, o_idx, Chem.BondType.DOUBLE)
            self._add_phenyl_ring(rwmol, p_idx)
        elif bridge == 'N-Ph':
            n_idx = rwmol.AddAtom(Chem.Atom(7))
            rwmol.AddBond(c1_idx, n_idx, Chem.BondType.SINGLE); rwmol.AddBond(c2_idx, n_idx, Chem.BondType.SINGLE)
            self._add_phenyl_ring(rwmol, n_idx)
        elif bridge == 'B-Ph':
            b_idx = rwmol.AddAtom(Chem.Atom(5))
            rwmol.AddBond(c1_idx, b_idx, Chem.BondType.SINGLE); rwmol.AddBond(c2_idx, b_idx, Chem.BondType.SINGLE)
            self._add_phenyl_ring(rwmol, b_idx)

        try:
            Chem.SanitizeMol(rwmol)
            return rwmol.GetMol(), True, "Rule2"
        except: return mol, False, "重排失败"

    def rule3_pi_extension(self, mol):
        rwmol = Chem.RWMol(mol)
        cands = []
        for b in rwmol.GetBonds():
            if b.GetIsAromatic():
                a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
                if a1.GetDegree() == 2 and a2.GetDegree() == 2 and a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                    cands.append(b)
        if not cands: return mol, False, "无边缘双键"
        
        bond = random.choice(cands)
        a1_idx, a2_idx = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        frag = random.choices(self.frag_types, weights=self.frag_weights, k=1)[0]
        
        if frag == 'Benzene': self._fuse_benzene(rwmol, a1_idx, a2_idx)
        elif frag == 'Furan': self._fuse_furan(rwmol, a1_idx, a2_idx)
        elif frag == 'Benzofuran': self._fuse_benzofuran(rwmol, a1_idx, a2_idx)
        elif frag == 'N-Ph-Indole': self._fuse_nphenylindole(rwmol, a1_idx, a2_idx)
            
        try:
            Chem.SanitizeMol(rwmol)
            return rwmol.GetMol(), True, "Rule3"
        except: return mol, False, "重排失败"

    def rule4_direct_coupling(self, mol):
        candidates, target_heteros = [], [5, 7, 8, 16, 34, 15]
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() in target_heteros:
                ipsos = sorted([n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIsAromatic()], key=lambda x: x.GetIdx())
                if len(ipsos) >= 2:
                    for i in range(len(ipsos)):
                        for j in range(i+1, len(ipsos)):
                            ip1, ip2 = ipsos[i], ipsos[j]
                            
                            in_same_ring = False
                            for ring in mol.GetRingInfo().AtomRings():
                                if atom.GetIdx() in ring and ip1.GetIdx() in ring and ip2.GetIdx() in ring:
                                    if len(ring) <= 7: in_same_ring = True; break
                            if in_same_ring: continue

                            orthos1 = [n for n in ip1.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIsAromatic() and n.GetTotalNumHs() > 0 and n.GetIdx() != atom.GetIdx()]
                            orthos2 = [n for n in ip2.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIsAromatic() and n.GetTotalNumHs() > 0 and n.GetIdx() != atom.GetIdx()]

                            bay_found = False
                            for o1 in orthos1:
                                for o2 in orthos2:
                                    if bay_found: break
                                    path = Chem.rdmolops.GetShortestPath(mol, o1.GetIdx(), o2.GetIdx())
                                    if len(path) == 5 and atom.GetIdx() in path:
                                        candidates.append((o1.GetIdx(), o2.GetIdx()))
                                        bay_found = True
        if not candidates: return mol, False, "无海湾"

        # 在演化循环中，一次只随机闭合一个海湾，保证演化的梯度和平滑性
        c1_idx, c2_idx = random.choice(candidates)
        rwmol = Chem.RWMol(mol)
        rwmol.AddBond(c1_idx, c2_idx, Chem.BondType.SINGLE)
        try:
            Chem.SanitizeMol(rwmol)
            return rwmol.GetMol(), True, "Rule4"
        except: return mol, False, "张力过大"

# ==========================================
# 自动化演化总控循环
# ==========================================
def run_dft_ml_generator(seed_smiles_list, num_iterations=2000, max_heavy_atoms=55):
    engine = MRMutationEngine()
    pool, unique_smiles, valid_seeds = [], set(), set()
    
    print("\n[1/4] 📥 初始化 MR-TADF 种子池...")
    for smi in seed_smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            can_smi = Chem.MolToSmiles(mol)
            pool.append(mol)
            unique_smiles.add(can_smi)
            valid_seeds.add(can_smi)
            print(f"  [+] 成功加载裸骨架: {can_smi}")
            
    if not pool: raise ValueError("无合法种子！")

    print(f"\n[2/4] 🚀 启动定向拓扑演化 (目标迭代: {num_iterations} 次)...")
    start_time = time.time()
    successful_mutations = 0
    
    for i in range(num_iterations):
        # 实时进度条
        if (i+1) % max(1, (num_iterations // 10)) == 0:
            elapsed = time.time() - start_time
            print(f"  ⏳ 进度: {i+1}/{num_iterations} | 已发现全新骨架: {len(unique_smiles) - len(valid_seeds)} 个 | 耗时: {elapsed:.1f}s")
            
        parent = random.choice(pool)
        
        # 权重分配：扩环负责造底物，掺杂调能级，桥连和闭环锁定刚性
        rules = [engine.rule1_n_doping, engine.rule2_planarization, engine.rule3_pi_extension, engine.rule4_direct_coupling]
        chosen_rule = random.choices(rules, weights=[0.05, 0.38, 0.38, 0.19], k=1)[0]
        
        child_mol, success, msg = chosen_rule(parent)
        
        if success and child_mol.GetNumHeavyAtoms() <= max_heavy_atoms:
            child_smi = Chem.MolToSmiles(child_mol)
            if child_smi not in unique_smiles:
                unique_smiles.add(child_smi)
                pool.append(child_mol) # 优质子代加入繁衍池
                successful_mutations += 1

    generated_smiles = list(unique_smiles - valid_seeds)
    print(f"\n[3/4] 🎉 演化结束！总计产生 {successful_mutations} 次有效跃迁，最终收获 {len(generated_smiles)} 个全新 MR 刚性骨架。")
    return generated_smiles

# ==========================================
# 执行入口与数据落盘
# ==========================================
if __name__ == "__main__":
    # 填入挑选的最基础的裸骨架 (Naked Cores)
    # DABNA-1 是 MR-TADF 领域的经典起点，结构简单且功能性强，非常适合演化发散
    seeds = [
        "C1=CC(=CC=C1)N1C2C(B3C4=CC=CC=C4N(C4C=CC=CC=4)C4C=CC=C1C=43)=CC=CC=2"
    ]
    
    # 执行大循环 (可以把 2000 改成 10000 跑个过瘾)
    novel_cores = run_dft_ml_generator(seed_smiles_list=seeds, num_iterations=2000, max_heavy_atoms=55)
    
    if novel_cores:
        print("\n[4/4] 💾 数据落盘与可视化准备...")
        
        # 1. 保存为 CSV，这是你们 DFT-ML 工作流最完美的输入格式
        csv_name = "AI_Generated_MR_Cores.csv"
        with open(csv_name, mode='w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(["ID", "SMILES", "Heavy_Atoms"])
            for idx, sm in enumerate(novel_cores):
                mol = Chem.MolFromSmiles(sm)
                writer.writerow([f"MR_Gen_{idx+1}", sm, mol.GetNumHeavyAtoms()])
        print(f"  -> SMILES 数据已导出至: {csv_name}")
        
        # 2. 渲染 20 个随机分子图谱供肉眼审核
        img_dir = "generator_results"
        if not os.path.exists(img_dir): os.makedirs(img_dir)
            
        display_smiles = random.sample(novel_cores, min(20, len(novel_cores)))
        mols_to_draw = [Chem.MolFromSmiles(sm) for sm in display_smiles]
        
        img = Draw.MolsToGridImage(
            mols_to_draw, molsPerRow=5, subImgSize=(400, 400),
            legends=[f"Gen_{i+1}" for i in range(len(mols_to_draw))],
            returnPNG=False
        )
        img_path = os.path.join(img_dir, "Top20_Showcase.png")
        img.save(img_path)
        print(f"  -> 高清骨架图谱已导出至: {img_path}")
        print("\n🚀 All Systems Go! 你们的专属 MR-TADF 发现引擎已完美竣工！")