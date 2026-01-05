import os
import gzip
import glob
import numpy as np
from Bio import PDB
from Bio.SeqUtils import seq1
import warnings

# 忽略 PDB 解析过程中的非致命警告
warnings.simplefilter('ignore', PDB.PDBExceptions.PDBConstructionWarning)

# ================= 配置区域 =================
# 输入输出路径
input_directory = os.path.expanduser("/workspace/zqlin/zexun_workbench/PeSTo/data/all_biounits")
output_fasta = "/workspace/zqlin/zexun_workbench/all_sequences_filtered.fasta"

# 过滤标准
MIN_LENGTH = 48         # 最小链长
MAX_RESOLUTION = 3.5    # 分辨率阈值
MAX_CA_DISTANCE = 4.2   # 断链阈值
# ===========================================

def has_large_gaps(chain_residues):
    """检查链内部是否有大的断裂"""
    if len(chain_residues) < 2: return False
    
    # 提取所有存在的 CA 原子坐标
    ca_coords = []
    for res in chain_residues:
        if 'CA' in res:
            ca_coords.append(res['CA'].coord)
            
    # 如果 CA 原子太少，无法计算 Gap，视作不完整
    if len(ca_coords) < 2: return True 
        
    coords = np.array(ca_coords)
    # 计算相邻 CA 距离
    diffs = coords[1:] - coords[:-1]
    dists = np.sqrt(np.sum(diffs**2, axis=1))
    
    # 只要有一个间距超过阈值，视为断链
    return np.any(dists > MAX_CA_DISTANCE)

def main():
    if not os.path.exists(input_directory):
        print(f"错误: 目录 {input_directory} 不存在")
        return

    search_path = os.path.join(input_directory, "**", "*.gz")
    all_gz_files = glob.glob(search_path, recursive=True)
    
    print(f"找到 {len(all_gz_files)} 个文件，开始处理...")
    print(f"策略: 分辨率>{MAX_RESOLUTION}丢弃(无分辨率则保留), 长度<{MIN_LENGTH}丢弃, 含UNK丢弃, 含Gap丢弃")

    parser = PDB.PDBParser(QUIET=True)
    
    stats = {
        "saved_chains": 0,       # 成功保存的链数量
        "skipped_resolution": 0, # 因分辨率差被丢弃的文件数
        "skipped_short": 0,      # 因太短被丢弃的链数
        "skipped_unk": 0,        # 因含 UNK 被丢弃的链数
        "skipped_gaps": 0,       # 因断链被丢弃的链数
        "error": 0               # 读取错误
    }

    with open(output_fasta, "w") as out_handle:
        for i, filepath in enumerate(all_gz_files):
            # 每处理 1000 个文件打印一次进度
            if i % 1000 == 0:
                print(f"进度: {i}/{len(all_gz_files)} | 已存链: {stats['saved_chains']} | 错误: {stats['error']}", end='\r')

            pdb_id = os.path.basename(filepath).split('.')[0]
            
            try:
                with gzip.open(filepath, 'rt') as f:
                    structure = parser.get_structure(pdb_id, f)
                    
                    # --- 关键修改 1: 分辨率检查策略变更 ---
                    # 只有当分辨率 "存在" 且 "大于阈值" 时才丢弃。
                    # 如果是 None，则默认保留！
                    resolution = structure.header.get('resolution')
                    if resolution is not None and resolution > MAX_RESOLUTION:
                        stats["skipped_resolution"] += 1
                        continue 

                    for model in structure:
                        for chain in model:
                            # 获取标准氨基酸
                            residues = [r for r in chain if PDB.is_aa(r, standard=True)]
                            
                            # --- 过滤: 链长 ---
                            if len(residues) < MIN_LENGTH:
                                stats["skipped_short"] += 1
                                continue
                            
                            # --- 过滤: UNK ---
                            if any(r.get_resname() == 'UNK' for r in residues):
                                stats["skipped_unk"] += 1
                                continue

                            # --- 过滤: 断链 (Gaps) ---
                            if has_large_gaps(residues):
                                stats["skipped_gaps"] += 1
                                continue

                            # --- 全部通过，保存序列 ---
                            seq_chars = []
                            for r in residues:
                                try:
                                    # 将三字母代码转单字母
                                    seq_chars.append(seq1(r.get_resname()))
                                except:
                                    seq_chars.append('X')
                            
                            sequence = "".join(seq_chars)
                            header = f">{pdb_id}_{chain.id}"
                            out_handle.write(f"{header}\n{sequence}\n")
                            
                            stats["saved_chains"] += 1
                            
            except Exception as e:
                # 某些文件可能损坏，跳过并记录
                stats["error"] += 1
                continue

    print(f"\n\n=== 处理完成 ===")
    print(f"总文件数: {len(all_gz_files)}")
    print(f"成功提取序列数: {stats['saved_chains']}")
    print(f"因分辨率差丢弃文件: {stats['skipped_resolution']}")
    print(f"丢弃短链 (<{MIN_LENGTH}): {stats['skipped_short']}")
    print(f"丢弃含 UNK 链: {stats['skipped_unk']}")
    print(f"丢弃断链 (Gaps): {stats['skipped_gaps']}")
    print(f"读取错误文件: {stats['error']}")
    print(f"结果已保存至: {output_fasta}")

if __name__ == "__main__":
    main()
