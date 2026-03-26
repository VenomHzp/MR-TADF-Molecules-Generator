from rdkit import Chem
import random

class Rule3Mutator:
    def __init__(self):
        # 法则 3：候选共轭片段及权重
        self.fragments = ['苯环 (Benzene)', '呋喃 (Furan)', '苯并呋喃 (Benzofuran)', 'N-苯基吲哚 (N-Ph-Indole)']
        self.weights = [0.55, 0.15, 0.15, 0.15]

    # ===================================================
    # 以下 4 个方法是底层的“外科缝合手术”，负责将片段完美接入骨架
    # ===================================================
    def _fuse_benzene(self, rwmol, a1_idx, a2_idx):
        """缝合片段 1：标准苯环扩增"""
        new_atoms = [rwmol.AddAtom(Chem.Atom(6)) for _ in range(4)]
        for idx in new_atoms: rwmol.GetAtomWithIdx(idx).SetIsAromatic(True)
        
        rwmol.AddBond(a1_idx, new_atoms[0], Chem.BondType.AROMATIC)
        rwmol.AddBond(new_atoms[0], new_atoms[1], Chem.BondType.AROMATIC)
        rwmol.AddBond(new_atoms[1], new_atoms[2], Chem.BondType.AROMATIC)
        rwmol.AddBond(new_atoms[2], new_atoms[3], Chem.BondType.AROMATIC)
        rwmol.AddBond(new_atoms[3], a2_idx, Chem.BondType.AROMATIC)

    def _fuse_furan(self, rwmol, a1_idx, a2_idx):
        """缝合片段 2：呋喃环 (O1*C=C*C=C1)"""
        o_idx = rwmol.AddAtom(Chem.Atom(8))
        c1_idx = rwmol.AddAtom(Chem.Atom(6))
        c2_idx = rwmol.AddAtom(Chem.Atom(6))
        for idx in [o_idx, c1_idx, c2_idx]: rwmol.GetAtomWithIdx(idx).SetIsAromatic(True)
        
        # 随机决定 O 原子的朝向
        if random.choice([True, False]):
            rwmol.AddBond(a1_idx, o_idx, Chem.BondType.AROMATIC)
            rwmol.AddBond(o_idx, c1_idx, Chem.BondType.AROMATIC)
            rwmol.AddBond(c1_idx, c2_idx, Chem.BondType.AROMATIC)
            rwmol.AddBond(c2_idx, a2_idx, Chem.BondType.AROMATIC)
        else:
            rwmol.AddBond(a2_idx, o_idx, Chem.BondType.AROMATIC)
            rwmol.AddBond(o_idx, c1_idx, Chem.BondType.AROMATIC)
            rwmol.AddBond(c1_idx, c2_idx, Chem.BondType.AROMATIC)
            rwmol.AddBond(c2_idx, a1_idx, Chem.BondType.AROMATIC)

    def _fuse_benzofuran(self, rwmol, a1_idx, a2_idx):
        """缝合片段 3：苯并呋喃"""
        o_idx = rwmol.AddAtom(Chem.Atom(8))
        rwmol.GetAtomWithIdx(o_idx).SetIsAromatic(True)
        
        bc = [rwmol.AddAtom(Chem.Atom(6)) for _ in range(6)]
        for idx in bc: rwmol.GetAtomWithIdx(idx).SetIsAromatic(True)
            
        for i in range(6):
            rwmol.AddBond(bc[i], bc[(i+1)%6], Chem.BondType.AROMATIC)
            
        if random.choice([True, False]):
            rwmol.AddBond(a1_idx, o_idx, Chem.BondType.AROMATIC)
            rwmol.AddBond(o_idx, bc[0], Chem.BondType.AROMATIC)
            rwmol.AddBond(bc[1], a2_idx, Chem.BondType.AROMATIC)
        else:
            rwmol.AddBond(a2_idx, o_idx, Chem.BondType.AROMATIC)
            rwmol.AddBond(o_idx, bc[0], Chem.BondType.AROMATIC)
            rwmol.AddBond(bc[1], a1_idx, Chem.BondType.AROMATIC)

    def _fuse_nphenylindole(self, rwmol, a1_idx, a2_idx):
        """缝合片段 4：N-苯基吲哚"""
        n_idx = rwmol.AddAtom(Chem.Atom(7))
        rwmol.GetAtomWithIdx(n_idx).SetIsAromatic(True)
        
        bc = [rwmol.AddAtom(Chem.Atom(6)) for _ in range(6)]
        ph = [rwmol.AddAtom(Chem.Atom(6)) for _ in range(6)]
        for idx in bc + ph: rwmol.GetAtomWithIdx(idx).SetIsAromatic(True)
            
        for i in range(6):
            rwmol.AddBond(bc[i], bc[(i+1)%6], Chem.BondType.AROMATIC)
            rwmol.AddBond(ph[i], ph[(i+1)%6], Chem.BondType.AROMATIC)
            
        # N-Ph 之间采用单键连接，允许空间旋转
        rwmol.AddBond(n_idx, ph[0], Chem.BondType.SINGLE)
        
        if random.choice([True, False]):
            rwmol.AddBond(a1_idx, n_idx, Chem.BondType.AROMATIC)
            rwmol.AddBond(n_idx, bc[0], Chem.BondType.AROMATIC)
            rwmol.AddBond(bc[1], a2_idx, Chem.BondType.AROMATIC)
        else:
            rwmol.AddBond(a2_idx, n_idx, Chem.BondType.AROMATIC)
            rwmol.AddBond(n_idx, bc[0], Chem.BondType.AROMATIC)
            rwmol.AddBond(bc[1], a1_idx, Chem.BondType.AROMATIC)


    # ===================================================
    # 法则 3 主逻辑：寻边与调用
    # ===================================================
    def rule3_pi_extension(self, mol):
        rwmol = Chem.RWMol(mol)
        
        # 识别外围空位：寻找连接的2个重原子数均为2的芳香碳-碳键
        candidate_bonds = []
        for b in rwmol.GetBonds():
            if b.GetIsAromatic():
                a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
                if a1.GetDegree() == 2 and a2.GetDegree() == 2 and a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                    candidate_bonds.append(b)
                    
        if not candidate_bonds:
            return None, "未找到可供扩增的边缘空位"
            
        bond = random.choice(candidate_bonds)
        a1_idx, a2_idx = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        
        # 按照权重抽取待插入的片段
        selected_fragment = random.choices(self.fragments, weights=self.weights, k=1)[0]
        
        # 执行对应的外科缝合手术
        if selected_fragment == '苯环 (Benzene)':
            self._fuse_benzene(rwmol, a1_idx, a2_idx)
        elif selected_fragment == '呋喃 (Furan)':
            self._fuse_furan(rwmol, a1_idx, a2_idx)
        elif selected_fragment == '苯并呋喃 (Benzofuran)':
            self._fuse_benzofuran(rwmol, a1_idx, a2_idx)
        elif selected_fragment == 'N-苯基吲哚 (N-Ph-Indole)':
            self._fuse_nphenylindole(rwmol, a1_idx, a2_idx)
            
        try:
            Chem.SanitizeMol(rwmol)
            return rwmol.GetMol(), selected_fragment
        except Exception as e:
            return None, selected_fragment


# ==========================================
# 运行测试
# ==========================================
if __name__ == "__main__":
    mutator = Rule3Mutator()
    
    # 经典的 DABNA-1 核心骨架
    dabna_smiles = "C1=CC(=CC=C1)N1C2C(B3C4=CC=CC=C4N(C4C=CC=CC=4)C4C=CC=C1C=43)=CC=CC=2"
    mol = Chem.MolFromSmiles(dabna_smiles)
    
    print(f"【初始状态】加载经典 DABNA 骨架: {dabna_smiles}")
    print("-" * 60)
    
    # 运行 8 次测试，观察概率分布和各种基团的插入效果
    for i in range(8):
        new_mol, fragment_used = mutator.rule3_pi_extension(mol)
        if new_mol:
            new_smiles = Chem.MolToSmiles(new_mol)
            print(f"[测试 {i+1}] 🎯 成功捕获边缘空位！")
            print(f"   拼装片段: {fragment_used}")
            print(f"   生成新型共轭扩增骨架: {new_smiles}\n")
        else:
            print(f"[测试 {i+1}] ❌ 扩增失败，尝试插入 {fragment_used} 时化合价未通过\n")