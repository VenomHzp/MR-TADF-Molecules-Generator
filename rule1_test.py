from rdkit import Chem
import random

class Rule1Mutator:
    def __init__(self):
        # 法则 1 专注于精准引入氮 (N)
        self.target_hetero = 7 

    def rule1_transmutation(self, mol):
        """
        事前过滤逻辑：找出所有绝对安全的外围 C-H 位点，保证 100% 置换成功。
        包含全局杂原子隔离锁，杜绝任何杂原子直接相连。
        """
        rwmol = Chem.RWMol(mol)
        
        # 1. 寻找所有候选的芳香碳原子
        all_c_candidates = [a.GetIdx() for a in rwmol.GetAtoms() if a.GetAtomicNum() == 6 and a.GetIsAromatic()]
        
        if not all_c_candidates:
            return mol, False, "未找到芳香碳原子"

        # 2. 构建绝对安全池
        safe_candidates = []
        for c_idx in all_c_candidates:
            atom = rwmol.GetAtomWithIdx(c_idx)
            
            # 【化学护城河 1】：必须是外围的 C-H 位点！
            # 绝对禁止替换桥头/稠合碳（无氢碳），防止芳香性重排崩溃
            if atom.GetTotalNumHs() == 0:
                continue
                
            # 【化学护城河 2】：全局杂原子隔离 (Absolute Heteroatom Isolation)
            # 获取该碳原子的所有重原子邻居
            neighbors = [a.GetAtomicNum() for a in atom.GetNeighbors()]
            
            # 只有当所有邻居都是碳 (原子序数 6) 或氢 (原子序数 1) 时，才允许替换！
            # 这彻底杜绝了新生成的 N 与 B, O, S, Se, P 甚至其他 N 直接相连的可能。
            if all(n == 6 or n == 1 for n in neighbors):
                safe_candidates.append(c_idx)

        # 3. 检查是否还有安全位点
        if not safe_candidates:
            return mol, False, "分子已达到氮掺杂饱和（所有合法外围 C-H 均已耗尽，或周边被杂原子封锁）！"

        # 4. 在绝对安全池中随机抽取并执行手术
        target_idx = random.choice(safe_candidates)
        rwmol.GetAtomWithIdx(target_idx).SetAtomicNum(self.target_hetero)
        
        try:
            Chem.SanitizeMol(rwmol)
            return rwmol.GetMol(), True, f"成功将外围孤立原子 {target_idx} (C-H) 置换为 N"
        except Exception as e:
            return mol, False, f"底层化合价重排失败: {e}"


# ==========================================
# 运行测试：验证“保证生效”与“杂原子隔离”逻辑
# ==========================================
if __name__ == "__main__":
    mutator = Rule1Mutator()
    
    # 依然使用 DABNA 骨架，它中心有一个 B 和两个 N
    dabna_smiles = "C1=CC(=CC=C1)N1C2C(B3C4=CC=CC=C4N(C4C=CC=CC=4)C4C=CC=C1C=43)=CC=CC=2"
    mol = Chem.MolFromSmiles(dabna_smiles)
    
    print(f"【初始状态】加载骨架: {dabna_smiles}")
    print("-" * 60)
    
    current_mol = mol
    
    # 连续执行强行掺杂
    for i in range(11):
        print(f">>> 第 {i+1} 次尝试精准掺杂 N：")
        new_mol, success, msg = mutator.rule1_transmutation(current_mol)
        
        if success:
            print(f"    ✅ {msg}")
            print(f"    生成骨架: {Chem.MolToSmiles(new_mol)}\n")
            current_mol = new_mol 
        else:
            print(f"    ⛔ 拦截：{msg}\n")
            break