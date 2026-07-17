#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
"""
File name    :   PanVariant
Created      :   2025/11/21
Author       :   Zeng Xiaojie 
Version      :   1.1
Contact      :   zengxiaojie@genomics.cn
License      :   Copyright (c) 2026, MGI.
"""

import os
from posixpath import isabs
import re
import sys
import logging
import argparse
import subprocess
from datetime import datetime
import signal

def add_path(path,path_all):
    path = path.rstrip('/')
    if os.path.isfile(path):
        path_dir = os.path.dirname(path)
        if path_dir not in path_all:
            path_all.add(path_dir)
    else:
        if path not in path_all:
            path_all.add(path)

def find_dir_symlinks(path,path_all):
    add_path(path,path_all)
    real_path = os.path.realpath(path)
    add_path(real_path,path_all)
    components = path.split(os.sep)
    for i in range(1, (len(components)+1)):
        sub_path = os.sep.join(components[:i])
        if os.path.islink(sub_path):
            resolved = os.readlink(sub_path)
            if not os.path.isabs(resolved):
                sub_path_up = os.path.dirname(sub_path)
                resolved = os.path.join(sub_path_up,resolved)
            add_path(resolved,path_all)
            find_dir_symlinks(resolved,path_all)

def find_file_symlinks(file_path,path_all):
    sampleList_dir=os.path.dirname(file_path)
    add_path(sampleList_dir,path_all)
    while os.path.islink(file_path):
        # 将符号链接解析为绝对路径
        resolved_path = os.path.abspath(os.readlink(file_path))
        resolved_dir = os.path.dirname(resolved_path)
        add_path(resolved_dir,path_all)
        # 获取符号链接所在的目录
        symlink_dir = os.path.dirname(file_path)
        # 添加到列表中
        add_path(symlink_dir,path_all)
        file_path = resolved_path
        # 如果解析后的路径是符号链接，继续循环
        if os.path.islink(file_path):
            continue
        # 否则，检查剩余的路径部分是否包含符号链接
        else:
            parent_dir = os.path.dirname(file_path)
            while parent_dir != symlink_dir:
                if os.path.islink(parent_dir):
                    add_path(parent_dir,path_all)
                    parent_dir = os.path.dirname(os.readlink(parent_dir))
                else:
                    break

def get_fq_path(samplelist,path_all):
    try:
        with open(samplelist, 'r', encoding='utf-8') as file:
            for line_num, line in enumerate(file, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                columns = line.split('\t')
                if len(columns) != 3:
                    print(f"Error: line {line_num} must have 3 tab-separated columns, got {len(columns)}.")
                    sys.exit(1)
                sample_name, read1_files, read2_files = columns
                
                read1_list = [path.strip() for path in read1_files.split(';') if path.strip()]
                read2_list = [path.strip() for path in read2_files.split(';') if path.strip()]
                if len(read1_list) != len(read2_list):
                    print(f"Error: line {line_num} has {len(read1_list)} read1 files but {len(read2_list)} read2 files.")
                    sys.exit(1)
                for i, (read1_path, read2_path) in enumerate(zip(read1_list, read2_list)):
                    if not os.path.isfile(read1_path):
                        print(f"Error: line {line_num} read1 file does not exist: {read1_path}")
                        sys.exit(1)
                    else:
                        find_file_symlinks(read1_path,path_all)
                    if not os.path.isfile(read2_path):
                        print(f"Error: line {line_num} read2 file does not exist: {read2_path}")
                        sys.exit(1)
                    else:
                        find_file_symlinks(read2_path,path_all)
                
                print(f"Line {line_num} validated: sample_name='{sample_name}', pairs={len(read1_list)}.")
    
    except FileNotFoundError:
        print(f"Error: sample list file does not exist: {samplelist}")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading sample list file: {e}")
        sys.exit(1)

def merge_paths(paths):
    abs_paths = [os.path.abspath(path) for path in paths]
    abs_paths.sort(key=lambda p: p.count(os.sep))
    merged_paths = []
    
    for path in abs_paths:
        if not any(os.path.commonpath([path, merged_path]) == merged_path for merged_path in merged_paths):
            merged_paths.append(path)
    return merged_paths

def main():
    LOG_FORMAT = "%(asctime)s %(levelname)s %(message)s "
    DATE_FORMAT = '%Y-%m-%d  %H:%M:%S %a ' 
    logging.basicConfig(level=logging.DEBUG,
                        format=LOG_FORMAT,
                        datefmt = DATE_FORMAT)
    script_path = os.path.dirname(os.path.realpath(__file__))
    panvariant_dir = '/usr/local/app/PanVariants'
    database_default = script_path + '/database'
    script_default = panvariant_dir + '/scripts'
    container_default = script_path + '/sifs/PanVariants.sif'
    nf_config_default = script_path + '/PanVariants.config'
    argparser = argparse.ArgumentParser(
        description='Run PanVariants pipeline by singularity container',
		prog='PanVariants',
		usage='./PanVariants [OPTIONS]',
        formatter_class=argparse.RawTextHelpFormatter)
    argparser.add_argument('-s','--samplelist', type=str, required=True,
    help="""Sample info, one or more lines with 3 columns split by tab:
        sample_name\tread1 path\tread2 path""")
    argparser.add_argument('-ex','--executor', type=str,default='blc',
    help="""The executor options, [blc|local].
        blc: Running pipeline using a Sun Grid Engine cluster
        local: Running pipeline using local machine""")
    argparser.add_argument('-o','--output',type=str,default='result',
                           help='Output directory, default [result]')
    argparser.add_argument('-ref','--reference',type=str,default='GRCh38',
                           help='The reference,[GRCh38 | CHM13], default [GRCh38]')
    argparser.add_argument('-rc','--run_cnv',type=str,default='yes',
                           help='Run CNV detection,[yes | no], default [yes]')
    argparser.add_argument('-rs','--run_sv',type=str,default='yes',
                           help='Run SV detection,[yes | no], default [yes]')
    argparser.add_argument('-rstr','--run_str',type=str,default='yes',
                           help='Run STR detection,[yes | no], default [yes]')
    argparser.add_argument('-md','--mark_dup',type=str,default='no',
                           help='Mark duplication,[yes | no], default [no]')
    argparser.add_argument('-ft','--fq_filter',type=str,default='no',
                           help='Run fq filter,[yes | no], default [no]')
    argparser.add_argument('-sb','--split_bam',type=str,default='no',
                           help='Split the bam file by chromosome and then run deepvariant by chromosome,[yes | no], default [no]')
    argparser.add_argument('-ad1','--adapter1',type=str,default='AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA',
                           help='The adapter1 sequences ,default [AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA]')
    argparser.add_argument('-ad2','--adapter2',type=str,default='AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG',
                           help='The adapter2 sequences ,default [AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG]')
    argparser.add_argument('-ra','--re_alignment',type=str,default='no',
                           help='Run re-alignment,[yes | no], default [no]')
    argparser.add_argument('-c','--nf_config',type=str,
                           default=nf_config_default,
                           help='config file nextflow')
    argparser.add_argument('-d','--db', type=str,
                           default=database_default,
                           help='database path')
    argparser.add_argument('-sp','--script', type=str,
                           default=script_default,
                           help='script path')
    argparser.add_argument('-f','--sif', type=str,
                           default=container_default,
                           help='singularity container path')
    argparser.add_argument('-scr','--scratch_tmp', type=str, required=False,
                           help='The path for nextflow scratch directive')
    argparser.add_argument('-q','--queue', type=str, default='none',
                           help='The queue name of the qsub command -q parameter. Note: It only takes effect when the executor option is set to blc.')
    argparser.add_argument('-pj','--project', type=str, default='none',
                           help='The project name of the qsub command -P parameter. Note: It only takes effect when the executor option is set to blc.')
    argparser.add_argument('-an','--ansi_log', type=str,
                           default='false',
                           help='Nextflow use ansi-log or not, [true|false], default [false]')
    args = argparser.parse_args()
    samplelist = os.path.realpath(args.samplelist)
    output = os.path.realpath(args.output)
    ref = args.reference
    run_cnv = args.run_cnv
    run_sv = args.run_sv
    run_str = args.run_str
    fq_filter = args.fq_filter
    split_bam = args.split_bam
    re_alignment = args.re_alignment
    mark_dup = args.mark_dup
    os.system(f'mkdir -p {output}/upload')
    nf_config = os.path.realpath(args.nf_config)
    database = os.path.realpath(args.db)
    script = os.path.realpath(args.script)
    script_dir = script if script != script_default else ''  
    container = os.path.realpath(args.sif)
    executor = args.executor
    queue = args.queue
    project = args.project
    ansi_log = args.ansi_log
    queue_para = f'"-q {queue}"' if queue != 'none' else '""'
    project_para = f'"-P {project}"' if project != 'none' else '""'
    current_time = datetime.now()
    formatted_time = current_time.strftime("%Y-%m-%d_%H%M%S")
    cwd_path = os.getcwd()
    scratch_tmp = args.scratch_tmp if args.scratch_tmp else ''
    scratch_test = scratch_tmp if args.scratch_tmp else 'false'
    nf_bind_list = [database,script_dir,output,cwd_path,script_path,scratch_tmp]
    all_nf_bind_list = set()
    get_fq_path(samplelist,all_nf_bind_list)
    for nf_b in nf_bind_list:
        find_dir_symlinks(nf_b,all_nf_bind_list)
    all_nf_bind_merge_list = merge_paths(all_nf_bind_list)
    all_nf_bind_str = ','.join(all_nf_bind_merge_list)
    runit = script_path + '/runit'
    # 配置 Nextflow 运行参数
    nextflow_config =f"""
    process{{
        container = \'{container}\'
        executor = \'{executor}\'
        errorStrategy = 'terminate'
        scratch = \'{scratch_test}\'
        stageOutMode = \'move\'
        maxRetries = 3
        maxErrors = 100
        beforeScript = 'echo "===vvv---HOST_INFO---vvv==="; uname -a; lscpu|grep -P "MHz|Model name:|^CPU:"; uptime; free -h; echo "===^^^---HOST_INFO---^^^==="; export MBP_USAGE_LOGDIR={cwd_path}/usage'
        shell = '{runit} -s -l `echo $(basename $(dirname $PWD))_$(basename $PWD)` -O 15 /bin/bash -euo pipefail'
    }}

    singularity{{
        enabled = true
        autoMounts = false
        singularity.runOptions = \"--env MBP_USAGE_LOGDIR={cwd_path}/usage,HOSTNAME=\\$HOSTNAME --no-home -C --pwd `pwd` -B `pwd -P` -B {all_nf_bind_str} -B {database}:/usr/local/app/PanVariants/Database -B `pwd -P`:/tmp\"
    }}
    """
    with open('nextflow.config','w') as nf:
        nf.write(nextflow_config)
    main_nf = script_path + '/main.nf'
    cmd = f"""
    mkdir -p tmp
    mkdir -p usage
    mkdir -p {output}
    rm -f trace.txt report.html timeline.html
    export DB={database}
    export SCRIPTS={script}
    export NXF_HOME={script_path}/NXF_HOME
    export NXF_OFFLINE=true
    export TMPDIR={cwd_path}/tmp
    {script_path}/nextflow run {main_nf} \
        -with-report report.html \
        -with-timeline timeline.html \
        -with-trace trace.txt \
        -ansi-log {ansi_log} \
        -profile {executor} \
        -resume \
        -c {nf_config} \
        --ref {ref} \
        --run_cnv {run_cnv} \
        --run_sv {run_sv} \
        --run_str {run_str} \
        --fq_filter {fq_filter} \
        --split_bam {split_bam} \
        --re_alignment {re_alignment} \
        --mark_dup {mark_dup} \
        --sampleList {samplelist} \
        --pipeline_database {database} \
        --runit {runit} \
        --queue {queue_para} \
        --project_para {project_para} \
        --outdir {output}
    """
    #print(cmd)
    signal.signal(signal.SIGTERM, handle_exit)
    signal.signal(signal.SIGINT, handle_exit)   # Ctrl+C

    # 启动子进程（需要创建进程组）
    global proc
    proc = subprocess.Popen(cmd,shell=True,start_new_session=True)
    return_code = proc.wait()
    if return_code != 0:
        sys.exit(return_code)

def handle_exit(signum, frame):
    # 终止子进程
    if 'proc' in globals() and proc.poll() is None:
        print("Received signal to terminate, killing subprocess...")
        os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
    exit(1)

if __name__ == '__main__' :
    main()

