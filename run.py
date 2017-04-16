# -*- coding: utf-8 -*-
import csv, sys, subprocess

def preprocessing(nodes, elems, title=''):
    num_node = len(nodes)
    num_elem = len(elems)

    alpha = 1
    beta = -1
    
    fid = open('input.dat', 'w')
    fid.write(title + '\n')
    fid.write('%d   %d   %d\n' % (num_elem, num_node, 0))
    fid.write('%10.6e  %10.6e\n' % (alpha, beta))
    fid.write('#Nodes:\n')
    for i, node in enumerate(nodes):
        fid.write(' %d     %20.12e     %20.12e  %20.12e  \n' % (i+1, node[0], node[1], node[2]))
    
    fid.write('#Elements:\n')
    for i, elem in enumerate(elems):
        fid.write(' %d     %d    %d     %d  %d  %20.12e \n' % (i+1, elem[0], elem[1], elem[2], elem[3], elem[4]))
        
    fid.close

def postprocessing():
     with open('output.dat', 'r') as output_file:
         data = [line.split() for line in output_file]
         output_data = {col[1]:list(col[2:]) for col in zip(*data)}
         # list to float
         Potential = [float(x) for x in output_data['Potential']]  
         Normal = [float(x) for x in output_data['Normal']]  
                      
         return Potential, Normal         

def read_csv(path):
    with open(path, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            yield [ float(i) for i in row ] 

if __name__ == '__main__':
    title = sys.argv[1] if len(sys.argv) > 1 else 'input'
    nodes = list(read_csv('nodes.csv'))
    elements = list(read_csv('elements.csv'))
    
    preprocessing(nodes, elements, title)
    
    #subprocess.call(['3D_Potential_CHBIE_FMM_64.exe'])

    #[phi, q] = postprocessing()