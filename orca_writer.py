import argparse
import os
import sys
import shutil



parser = argparse.ArgumentParser()

parser.add_argument('file', help='Scan log file(s)', nargs='+')



args = parser.parse_args()



INPUT = '''
! cam-b3lyp 6-311++g(d,p) CPCM(hexane)

% PAL
    NPROCS 44
end


%tddft
mode stddft
end

*xyz 0 1 
{atoms}
end
'''

def create_input(filename):

    fn = f'{filename.split(".")[0]}.inp'
    directory = filename.split('.')[0]
    with open(filename) as f:
        atoms = ''.join(f.readlines()[2:])
    with open(fn, 'w') as f:
        f.write(INPUT.format(atoms=atoms).strip()+'\n')
    
    os.mkdir(os.path.join('inputs_orca', directory))
    os.rename(filename, os.path.join('inputs_orca', directory, filename))
    os.rename(fn, os.path.join('inputs_orca', directory, fn))



if __name__=='__main__':

    os.mkdir('inputs_orca')
    for file in args.file:
        create_input(file)