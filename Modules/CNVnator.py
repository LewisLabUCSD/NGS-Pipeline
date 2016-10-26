import sarge,sys
def cnvnator(bam,out_root,otherParameters=['']):
    '''run cnvnator'''
    # 1. predict CNV regions
    cmd = ('cnvnator -root {out} -tree {bam} {other}').format(out=out_root,bam=bam,
                                                        other=' '.join(otherParameters))
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    # 2. generate histograme
    bin = '100'
    path = '/'.join(out_root.split('/')[:-1])
    cmd = ('cnvnator -root {root} -his {bin} -d {dir}').format(root=out_root,bin=bin,dir=path)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    # 3. calculate statistics
    cmd = ('cnvnaotr -root {root} -stat {bin}').format(root=out_root,bin=bin)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    # 4. rd signal partition
    cmd = ('cnvnator -root {root} -partition {bin}').format(root=out_root,bin=bin)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    # 5. call CNV
    out = 'cnv/cnv_res.txt'
    cmd = ('cnvnator -root {root} -call {bin} > {out}').format(root=out_root,bin=bin,out=out)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    
    
    