"""
this script is used to convert a .txt file containing gaits to
.ino file which is the arduino code that runs the robot

1. Read the .ino file
2. Generate place holders
3. Read .txt file
4. Generate the .ino file

Run it using the following template:
python genarduino.py  -i gait2018_3.ino -t forward.txt -o gait2018_forward.ino
"""

import argparse
import os, re, sys
from pprint import pprint

# var_map = {
# 	'F,L,_,a,n,g,_,l,a,t,e,r,a,l': 'FL_ang_lateral',
# 	'F,L,_,a,n,g,_,v,e,r,t': 'FL_ang_vert',
# 	'F,R,_,a,n,g,_,l,a,t,e,r,a,l': 'FR_ang_lateral',
# 	'F,R,_,a,n,g,_,v,e,r,t': 'FR_ang_vert',
# 	'H,L,_,a,n,g,_,l,a,t,e,r,a,l': 'HL_ang_lateral',
# 	'H,L,_,a,n,g,_,v,e,r,t': 'HL_ang_vert',
# 	'H,R,_,a,n,g,_,l,a,t,e,r,a,l': 'HR_ang_lateral',
# 	'H,R,_,a,n,g,_,v,e,r,t': 'HR_ang_vert',
# 	'b,o,d,y,_,l,a,t,e,r,a,l,_,1': 'body_lateral1',
# 	'b,o,d,y,_,l,a,t,e,r,a,l,_,2': 'body_lateral2',
# 	'b,o,d,y,_,v,e,r,t': 'body_vert'
# }

bad_var = ['body_lateral_1', 'body_lateral_2']

def read_arguments():
	parser = argparse.ArgumentParser(description='Generate the ino file with correct gaits')
	parser.add_argument('-i','--inp_ino', help='Input ino file', required=True)
	parser.add_argument('-t','--inp_txt', help='Input txt file', required=True)
	parser.add_argument('-o','--out_ino', help='Output txt file', required=True)
	args = vars(parser.parse_args())

	inp_ino = args['inp_ino']
	inp_txt = args['inp_txt']
	out_ino = args['out_ino']

	return inp_ino, inp_txt, out_ino

def extract_var(line):
	result = re.search('int (.*)\[', line)
	return result.group(1)

def read_ino_file(input_filename):
	saved = []
	temp_str = ''
	pattern = re.compile("int [a-zA-Z_0-9]*\[")
	with open(input_filename) as file:
		for line in file:
			if pattern.match(line):
				saved.append(temp_str)
				temp_str = ''
				var = extract_var(line)
				saved.append(var)
				continue
			temp_str += line
	saved.append(temp_str)
	return saved

def read_txt_file(input_filename):
	flag = 'var'
	variable = ''
	place_holder = {}
	with open(input_filename) as file:
		for line in file:
			if flag == 'var':
				variable = line[:-1].replace(',','')
				if variable in bad_var:
					variable = variable[:-2] + variable[-1:]
				flag = 'val'
			elif flag == 'val':
				value = line[:-1]
				place_holder[variable] = value
				flag = 'var'
	return place_holder

def format_line(var, val):
	string = 'int '+var+'[241] = {'+val+'};'
	return string

def gen_save_output(output_filename, content, place_holder):
	with open(output_filename, 'w') as file:
		for line in content:
			if line in place_holder:
				print (format_line(line, place_holder[line]))
			else:
				print (line)

if __name__ == '__main__':
	input_ino, input_txt, output_ino = read_arguments()
	ino_content = read_ino_file(input_ino)
	place_holder = read_txt_file(input_txt)
	output_content = gen_save_output(output_ino, ino_content, place_holder)
	# pprint (place_holder)