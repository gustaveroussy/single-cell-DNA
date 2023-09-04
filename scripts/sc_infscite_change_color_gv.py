import re
import random
import seaborn as sns
import argparse

parser = argparse.ArgumentParser(description='change color of graphviz tree')

parser.add_argument("--input_gv", help="input gv path file",required=True,action='store')
parser.add_argument("--output_gv", help="output gv path file",required=True,action='store')

args = parser.parse_args()

#def generate_random_color():
#    return '#{:06x}'.format(random.randint(0x00d0FF,0x00FF00))

with open(args.input_gv,"r") as stream:
    list_line=[i_line for i_line in stream]
    
pattern= "^\d*\[.*"

res = [x for x in list_line if re.search(pattern, x)]
res=[i.split('"') for i in res]
#color=[generate_random_color() for i in range(len(res))]

color=sns.color_palette("pastel",len(res))
color=list(color.as_hex())

for i in range(0,len(res)):
    res[i][2]=',color="'+color[i]+'"];\n'
    
res=[res[i][0]+'"'+res[i][1]+'"'+res[i][2] for i in range(0,len(res))]

for i in range(0,len(res)):
    list_line[i+2]=res[i]

res_2=list_line[len(res)+2:(len(res)*2)+1]
res_2=[i.split(";") for i in res_2]
res_2=[res_2[i][0]+" [dir=none style=dashed weight=1 penwidth=5];"+res_2[i][1] for i in range(0,len(res_2))]

reverse_line_list=list_line[::-1]

for i in range(1,len(res_2)+1):
    reverse_line_list[i]=res_2[i-1]
    
list_line=reverse_line_list[::-1]
    
with open(args.output_gv,"w") as stream:
    for item in list_line:
        stream.write(item)