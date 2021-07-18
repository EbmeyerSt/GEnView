import os
import csv
import math
from types import new_class

path = 'C:/Users/xcoero/Documents/Projects/scripts/GEnView/'
os.chdir(path)


def make_text(final_tree):
    if os.path.exists("index_tree.txt"):
        os.remove("index_tree.txt")
    f = open("index_tree.txt", "a")
    for row in final_tree:
        f.write(str(row))
        f.write('\n')
    f.close()

def get_parent(child):
    t = list(child)
    f = 0
    for v in reversed(t):
        if v.isnumeric():
            f += 1
        else:
            f += 1
            break
    k = ''.join(t)
    k = k[:-f]
    return k

def get_last_letter(string):
    t = list(string)
    res = []
    f = 0
    for v in reversed(t):
        if v.isnumeric():
            res.insert(0, v)
        else:
            break
    k = ''.join(res)
    return k

def create_tree_index(tree_file_path):
    with open(tree_file_path) as tree:
        completed = False
        array = []
        final_tree = []
        #Extract first line of tree file which contains the tree information
        for line in tree: 
            #Remove all bootstrap values from data      
            i = 0
            line = list(line)
            while True:
                if line[i] == ')' and line[i+6] == ':':
                    del line[i+1:i+6]
                i += 1
                if i == len(line)-6:
                    break 
            line = ''.join(line)
            #print(line)
            used_ids = []
            array.extend(line)
            current_id = '1'
            #Iterate through each character in tree file line
            i = 0
            while i < len(array):
                #"(" indicates a new branch and therefore creates a branch array
                if array[i] == ')' and array[i+1] == ';':
                    completed = True
                    break
                elif array[i] == "(":   
                    #Adjust id's and add to ID list
                    id_ext = 1
                    while True:
                        tmp_id = current_id + '.' + str(id_ext)
                        if tmp_id in used_ids:
                            id_ext += 1
                        else:
                            current_id = tmp_id
                            break    

                    if len(used_ids) == 0:
                        final_tree.append([current_id])
                    else:
                        final_tree.append([current_id]) 
                    used_ids.append(current_id) 
                #Skip bootsrtap values
                elif array[i] == ')' and array[i+6] == ':':
                    i += 5
                elif array[i] == ":" and array[i-1] != ")":
                    #Adjust id's and add to ID list
                    id_ext = 1
                    while True:
                        tmp_id = current_id + '.' + str(id_ext)
                        if tmp_id in used_ids:
                            id_ext += 1
                        else:
                            break
                    final_tree.append([tmp_id])
                    used_ids.append(tmp_id)    
                #Marks end of branch and all sub branches so length of this branch must be added
                elif array[i] == ":" and array[i-1] == ")":              
                    current_id = get_parent(current_id)
                i += 1
            if completed == True:
                break

    result = []
    for item in final_tree:
        l = len(item[0])/2+0.5
        result.append(l)
    res = max(result)
    #Organize tree branches
    p = 1
    new_tree = []
    while p <= res:
        for item in final_tree:
            if int((len(item[0])/2+0.5)) == p:
                parent = get_parent(item[0])
                try:
                    n = new_tree.index([parent])    
                except:
                    n = 0                            
                last_letter = int(get_last_letter(item[0]))
                if last_letter == 1:
                    new_tree.insert(n, item)
                elif last_letter == 2:
                    new_tree.insert(n+1, item)
                else:
                    last_sibling = parent + '.' + str(int(last_letter)-1)
                    n = new_tree.index([last_sibling]) 
                    new_tree.insert(n+1, item)
                    new_tree.insert(n+1, ['-'])
        p += 1            

    #Create list of orderd organism names
    with open(tree_file_path) as tree:
        final_order = []
        for line in tree: 
            line = list(line)
            i = 0
            while i < (len(line)-1):
                if line[i] != '(' and line[i-1] == '(':
                    name = ''
                    while True:
                        name += line[i]
                        i += 1
                        if line[i] == ':':
                            final_order.append(name)
                            break
                elif line[i-1] == ',' and line[i] != '(':
                    name = ''
                    while True:
                        name += line[i]
                        i += 1
                        if line[i] == ':':
                            final_order.append(name)
                            break
                i += 1




    space = [0]*5
    lne = [1]*5
    #Assign distances to branches and nodes
    for item in new_tree:
        #Check for if branch or node
        again = 0
        for branch in new_tree:
            if branch[0].startswith(item[0]):
                again += 1
        if item[0] == '-': 
            item.extend(space)
        elif item[0] == '1.1':
            item.extend(space)
        elif again <= 1:
            #If node
            a = int((len(item[0])/2+0.5))-2
            pre = space*a
            item.extend(pre)
            item.extend(lne)
        else:
            #if branch
            a = int((len(item[0])/2+0.5))-2
            pre = space*a
            item.extend(pre)
            item.extend(lne)

    #Add white spacer before and after lists
    new_tree.insert(0, [0]*len(new_tree[0]))
    new_tree[0][0] = '-'
    new_tree.append([0]*len(new_tree[0]))
    new_tree[len(new_tree)-1][0] = '-'

    #Adding spacers between list elements where required
    i = 0 
    while i < len(new_tree[:-1]):
        n = i + 1
        parent_i = get_parent(new_tree[i][0])
        parent_n = get_parent(new_tree[n][0])
        if parent_i == parent_n:
            spacer = [0]*(len(new_tree[i]))
            spacer[0] = '-'
            new_tree.insert(n, spacer)
        i += 1
    long = len(max(new_tree, key=len)) + 1
    z = 0
    while z < len(new_tree):
        if len(new_tree[z]) <= long:
            val = long - len(new_tree[z])
            p = 1
            while p <= val:
                new_tree[z].append(0)
                p += 1
        z += 1
    i = 0
    while i < len(new_tree):
        new_tree[i].insert(1, 0)
        i += 1
    i = 0
    while i < len(new_tree):
        new_tree[i].insert(1, 0)
        i += 1

    #Making lengths of list names for each branch and node the same length
    lenghts = []
    for row in new_tree:
        lenghts.append(len(str(row[0])))
    m = max(lenghts)
    i = 0
    while i < len(new_tree):
        add = m - len(str(new_tree[i][0]))
        if add > 0:
            new_tree[i][0] = str(new_tree[i][0]) + '-'*add
        i += 1


    #Center horizontal lines
    #Get max number of levels
    mx = 0
    for branch in new_tree[1:-1]:
        name = branch[0].strip('-')
        cur_level = name.count('.')
        if cur_level > mx:
            mx = cur_level

    shift_index = []
    while mx >= 1:
        print(mx)
        cur_row = 1
        for b in new_tree[1:-1]:
            name = b[0].strip('-')
            cur_level = name.count('.')
            if cur_level == mx:
                if len(name) > 0:
                    #Intercept here to remove and replace 1's
                    r_start = (cur_level-1)*5 + 3
                    r_end = r_start + 4
                    again = 0
                    #See if b is branch or node
                    pos_child = name + '.1'
                    for branch in new_tree[1:-1]:
                        if branch[0].startswith(pos_child):
                                again += 1
                    #IF branch get all direct children and count distance.
                    if again > 0:
                        ext = 1
                        tmp_row = 1
                        for branch in new_tree[1:-1]:
                            child = name + '.' + str(ext)
                            if branch[0].strip('-') == child:
                                if ext == 1:
                                    if child in shift_index:
                                        shift = shift_index.index(child)
                                        first_row = shift_index[shift+2]
                                    else:
                                        first_row = tmp_row                                        
                                    ext += 1
                                else:
                                    if child in shift_index:
                                        shift = shift_index.index(child)
                                        last_row = shift_index[shift+2]
                                    else:
                                        last_row = tmp_row
                                    ext += 1
                            tmp_row += 1 
                        #Remove existing 1's in range and replace in midpoint
                        rel_mid = math.ceil(((last_row-first_row)/2))
                        actual_mid_row = first_row+rel_mid                        
                        if cur_row != actual_mid_row:
                            shift_index.append(name)
                            shift_index.append(cur_row)
                            shift_index.append(actual_mid_row)
                            i = r_start
                            while i <= r_end:
                                new_tree[cur_row][i] = 0
                                i += 1                            
                            i = r_start
                            while i <= r_end:
                                new_tree[actual_mid_row][i] = 1
                                i += 1
                        #Add vertical lines:  
                        ver_line = r_end + 1
                        k = first_row
                        while k <= last_row:
                            new_tree[k][ver_line] = 1
                            k += 1
            cur_row += 1
        mx -= 1

    #Vertically connect horizontal lines
    # i = 1
    # while i < len(new_tree[0]):
    #     n = 0
    #     while n < len(new_tree[1:-1]):
    #         if new_tree[n][i] == 1 and new_tree[n][i-1] == 0:
    #             parent = new_tree[n][0].strip('-')
    #             parent = get_parent(parent)                
    #             cur = new_tree[n][0]
    #             cur = cur.strip('-')
    #             ll = int(get_last_letter(cur))
    #             number = 0
    #             while True:
    #                 ll += 1
    #                 tmp = parent + '.' + str(ll)
    #                 cont = False
    #                 for branch in new_tree[1:-1]:
    #                     if branch[0].strip('-') == tmp:
    #                         number += 1
    #                         cont = True
    #                 if cont == False:
    #                     break
    #             ll = int(get_last_letter(cur)) + 1
    #             while True:
    #                 if number == 0:
    #                     break
    #                 else:
    #                     tmp = parent + '.' + str(ll)
    #                     new_tree[n][i] = 1
    #                     if new_tree[n][0].strip('-') == tmp:
    #                         number -= 1
    #                         ll += 1
    #                     n += 1
    #         n += 1
    #     i += 1

            
    #Add stipple lines to nodes
    i = 1
    while i < len(new_tree[:-1]):
        n = len(new_tree[0])-3
        while n >= 0:
            if new_tree[i][n-2] == 1:
                break
            else:
                new_tree[i][n] = 1
                new_tree[i][n-1] = 1
            n -= 4
        i += 2


    make_text(new_tree)
    return new_tree

meta_data_path = "annotation_meta.csv"
tree_file_path = "iscr_contexts.unique.txt"

final_tree = create_tree_index(tree_file_path)

#Create string of phylogenetic tree HTML
tree_string = ''
x_sep = round(300/len(final_tree[0][1:]),2)


tree_string += '<svg class="tree l_tree" width="300" height="12000">'
#Get horizontal lines
y = 17
for row in final_tree[1:-1]:
    x1 = 0
    x2 = 0
    xval = 0
    i = 1
    while i < len(row):
        n = i + 1
        p = i - 1
        if row[i] == 1 and row[n] == 1 and row[p] == 0:
            x1 = xval
        elif row[i] == 1 and row[n] == 0 and row[p] == 1:
            x2 = xval
            tree_string += '<line class="hor" x1="'+ str(x1) +'" y1="'+ str(y) +'" x2="'+ str(x2) +'" y2="'+ str(y) +'" stroke="black" stroke-width="1"/>'
        xval += x_sep
        i += 1
    y += 16

#get vertical lines
i = 1
xval = 0
while i < len(final_tree[0]):
    y1 = 0
    y2 = 0
    y = 17
    k = 1    
    while k < (len(final_tree)-1):
        if final_tree[k][i] == 1 and final_tree[k+1][i] == 1 and final_tree[k-1][i] == 0:
            y1 = y
        elif final_tree[k][i] == 1 and final_tree[k+1][i] == 0 and final_tree[k-1][i] == 1:
            y2 = y
            if y1 > 0 and y2 > 0:
                tree_string += '<line class="ver" x1="'+ str(xval) +'" y1="'+ str(y1) +'" x2="'+ str(xval) +'" y2="'+ str(y2) +'"/>'
        k += 1    
        y += 16        
    xval += x_sep
    i += 1
tree_string += '</svg>'


tree_string += '<svg class="tree m_tree hidden" width="300" height="6000">'
    #Get horizontal lines
y = 10
for row in final_tree[1:-1]:
    x1 = 0
    x2 = 0
    xval = 0
    i = 1
    while i < len(row):
        n = i + 1
        p = i - 1
        if row[i] == 1 and row[n] == 1 and row[p] == 0:
            x1 = xval
        elif row[i] == 1 and row[n] == 0 and row[p] == 1:
            x2 = xval
            tree_string += '<line class="hor" x1="'+ str(x1) +'" y1="'+ str(y) +'" x2="'+ str(x2) +'" y2="'+ str(y) +'" stroke="black" stroke-width="1"/>'
        xval += x_sep
        i += 1
    y += 9
#get vertical lines
i = 1
xval = 0
while i < len(final_tree[0]):
    y1 = 0
    y2 = 0
    y = 10
    k = 1
    while k < (len(final_tree)-1):
        p = k - 1
        n = k + 1
        if final_tree[k][i] == 1 and final_tree[n][i] == 1 and final_tree[p][i] == 0:
            y1 = y
        elif final_tree[k][i] == 1 and final_tree[n][i] == 0 and final_tree[p][i] == 1:
            y2 = y
            if y1 > 0 and y2 > 0:
                tree_string += '<line class="ver" x1="'+ str(xval) +'" y1="'+ str(y1) +'" x2="'+ str(xval) +'" y2="'+ str(y2) +'"/>'
        k += 1    
        y += 9       
    xval += x_sep
    i += 1
tree_string += '</svg>'


tree_string += '<svg class="tree s_tree hidden" width="300" height="3000">'
    #Get horizontal lines
y = 5
for row in final_tree[1:-1]:
    x1 = 0
    x2 = 0
    xval = 0
    i = 1
    while i < len(row):
        n = i + 1
        p = i - 1
        if row[i] == 1 and row[n] == 1 and row[p] == 0:
            x1 = xval
        elif row[i] == 1 and row[n] == 0 and row[p] == 1:
            x2 = xval
            tree_string += '<line class="hor" x1="'+ str(x1) +'" y1="'+ str(y) +'" x2="'+ str(x2) +'" y2="'+ str(y) +'" stroke="black" stroke-width="1"/>'
        xval += x_sep
        i += 1
    y += 4
#get vertical lines
i = 1
xval = 0
while i < len(final_tree[0]):
    y1 = 0
    y2 = 0
    y = 5
    k = 1
    while k < (len(final_tree)-1):
        p = k - 1
        n = k + 1
        if final_tree[k][i] == 1 and final_tree[n][i] == 1 and final_tree[p][i] == 0:
            y1 = y
        elif final_tree[k][i] == 1 and final_tree[n][i] == 0 and final_tree[p][i] == 1:
            y2 = y
            if y1 > 0 and y2 > 0:
                tree_string += '<line class="ver" x1="'+ str(xval) +'" y1="'+ str(y1) +'" x2="'+ str(xval) +'" y2="'+ str(y2) +'"/>'
        k += 1    
        y += 4    
    xval += x_sep
    i += 1
tree_string += '</svg>'





with open(meta_data_path, newline = '') as org_file:
    reader = csv.reader(org_file, delimiter='\t')
    values = []
    for row in reader:
        if len(row) > 1:
            values.append(int(row[2]))


            
factor = 1000/max(values)


name = False
header = ""
start_next = ""
string = ""
run = 0
id = 0
x = True

tnps=['iscr', 'transpos', 'tnp', 'insertion']
ints=['inti', 'integrase', 'xerc', 'xerd']
mobiles=['secretion', 'mobiliza', 'moba', 'mobb', 'mobc', 'mobl', 'plasmid', 'relaxase', 'conjugation', 'type iv']
res=['lactam', 'aminoglyco', 'fluoroquinolo', 'tetracyclin', 'macrolid', 'carbapenem']

line_no = 2
line_no_2 = line_no + 2
gene_no_row_2 = 1
string += '<div class="sequences" style="grid-column-start: 2; grid-column-end: 3;grid-row-start: 1;grid-row-end: 2;">'


with open(meta_data_path, newline = '') as file:
    reader = csv.reader(file, delimiter='\t')
    for row in reader:
        print(row)

with open(meta_data_path, newline = '') as file:
    reader = csv.reader(file, delimiter='\t')
    for row in reader:
        if len(row) > 0:
            if len(row) == 1 and "__________" in row[0] and name == False: 
                name = True
                run += 1
            elif name == True:
                header = row[0]
                if x == True:
                    x = False
                else:
                    string += '</div>'
                string += '<div class="grid_container">'
                string += '<div id="' + str(id) + '_line" class="line" style="grid-column: 1 / 1000;grid-row: 2 / 4;background-color: rgba(0, 0, 0);opacity: 0.2;"></div>'
                line_no = line_no + 4
                line_no_2 = line_no + 2
                name = False
            else:
                name = False
            if len(row) == 1 and "---------------" in row[0]:
                gene_no_row = gene_no_row_2
                gene_no_row_2 = gene_no_row + 4
                start_next = True
            elif len(row) == 1 and "---------------" not in row[0]:
                start_next = False
            elif start_next == True:
                start = math.ceil(int(row[1])*factor)
                end = math.ceil(int(row[2])*factor)

                if any(keyword in row[0].lower() for keyword in tnps):
                    color='violet'
                    color = 'rgb(0, 181, 253)'
                elif any(keyword in row[0].lower() for keyword in ints):
                    color='yellow'
                    color = 'rgb(253, 228, 0)'
                elif any(keyword in row[0].lower() for keyword in mobiles):
                    color='green'
                    color = 'rgba(0, 255, 13)'
                elif any(keyword in row[0].lower() for keyword in res):
                    color='DodgerBlue'
                    color = 'rgba(0, 183, 255)'
                elif 'hypothetical' in row[0].lower():
                    color='grey'
                    color = 'rgb(153, 153, 153)'

                if row[3] == '+' or row[3] == '-':
                    if row[3] == '+':
                        string += '<div id="'+ str(id) +'_gene" class="R all" style="grid-column: '+ str(start) +' / '+ str(end) +';grid-row: 1 / 5;background-color: '+ color +';margin-right: 15px;opacity: 0.7;"><p></p></div>'
                        string += '<div class="AR right" id="'+ str(id) +'_gene_arrow" style="grid-column: '+ str(end) +' / '+ str(end) +';grid-row: 1 / 5;--my-color-var: '+ color +';opacity: 0.7;"></div>'
                    elif row[3] == '-':
                        string += '<div class="AL left" id="'+ str(id) +'_gene_arrow" style="grid-column: '+ str(start) +' / '+ str(start) +';grid-row: 1 / 5;--my-color-var: '+ color +';opacity: 0.7;"></div>'
                        string += '<div id="'+ str(id) +'_gene" class="L all" style="grid-column: '+ str(start) +' / '+ str(end) +';grid-row: 1 / 5;background-color: '+ color +';margin-left: 15px;opacity: 0.7;"><p></p></div>'

                    string += '<div id="'+ str(id) +'_gene_info" class="hidden info_box"><button class="exit">X</button><p><strong>Organism:</strong> Citrobacter farmeri</p><p><strong>Accession:</strong> H23409823.1</p><p><strong>Gene:</strong> '+ row[0] +'</p><p><strong>Start:</strong> '+ str(row[1]) +'</p><p><strong>End:</strong> '+ str(row[2]) +'</p><textarea id="'+ str(id) +'_gene_sequence">'+ header +'&#13;&#10;'+ row[4] +'</textarea><button class="copy btn" onclick="copy_to_clipboard("'+ str(id) +'_gene_sequence");">Copy</button></div>'           
            id += 1   
    string += '</div>'  

            



    

      
html = """
<html>
<head>
    <title>GenView Output</title>
    <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js"></script>
</head>

<style>

.space {
    position: relative;
    height: 50px;
}

.navigation {
    top: 0;
    border: none;
    outline: none;
    position: fixed;
    width: 100%;
    height: 50px;
    background-color: rgb(131, 131, 131);
    margin-bottom: 60px;
    display: inline_block;
    overflow: hidden;
    z-index: 100;
    box-shadow: 0px 5px 5px rgb(131, 131, 131, 0.5);
}

.active {
    background-color: rgb(70, 181, 201) !important; 
}

.navigation button {
    position: relative;
    width: 100px;
    background-color: rgb(131, 131, 131);
    display: inline-block;
    border: none;
    outline: none;
    transition: 0.3s;
    color: white;
    cursor: pointer; 
    height: 100%;
}

.navigation button:hover {
    background-color: rgb(70, 181, 201);
    transition: 0.3s;
}



.navigation button:active {
    background-color: rgb(131, 131, 131);
    transition: 0.3s;
}

line {
    stroke: black;
    stroke-width: 1;
}

.phylo_div {
    height: 100%;
}

.phylo_div div {
    display: inline-block;
    height: 10px;
    padding: 0px;
}

.large_grid_container {
    display: grid;
    height: 100%;
    grid-template-columns: 300px auto;
}

.grid_container {
  margin-top: 2px;
  display: grid;
  grid-row-gap: 2px;
  grid-template-columns: repeat(1000, 1fr);
}

.grid_container p {
    margin-top: 6px;
    margin-bottom: 0px;
    padding: none;
    cursor: pointer;  
}

.all {
    transition: 0.3s;
    text-align: center;
    margin: none;
    z-index: 2;
    transition: 0.3s;
    height: 30px;
}

.line {
    z-index: 1;
    transition: 0.3s;
    cursor: pointer; 
    height: 10px; 
}


.left {    
    position: relative;
    outline: none;
    width: -webkit-calc(100% - 15px);
    width:    -moz-calc(100% - 15px);
    width:         calc(100% - 15px);
    z-index: 2;
    transition: 0.3s;
}
.left:after {
    content: '';
    position: absolute;
    top: 0; 
    border-top: 15px solid transparent;
    border-right: 15px solid var(--my-color-var);
    border-bottom: 15px solid transparent;
    width: 0;
}


.left_s {    
    position: relative;
    outline: none;
    width: -webkit-calc(100% - 15px);
    width:    -moz-calc(100% - 15px);
    width:         calc(100% - 15px);
    z-index: 2;
    transition: 0.3s;
}
.left_s:after {
    content: '';
    position: absolute;
    top: 0; 
    border-top: 3px solid transparent;
    border-right: 6px solid var(--my-color-var);
    border-bottom: 3px solid transparent;
    width: 0;
}
.left_m {    
    position: relative;
    outline: none;
    width: -webkit-calc(100% - 15px);
    width:    -moz-calc(100% - 15px);
    width:         calc(100% - 15px);
    z-index: 2;
    transition: 0.3s;
}
.left_m:after {
    content: '';
    position: absolute;
    top: 0; 
    border-top: 8px solid transparent;
    border-right: 8px solid var(--my-color-var);
    border-bottom: 8px solid transparent;
    width: 0;
}


.right {    
    position: relative;
    outline: none;
    width: -webkit-calc(100% - 15px);
    width:    -moz-calc(100% - 15px);
    width:         calc(100% - 15px);
    left: -15px;
    z-index: 2;
    transition: 0.3s;
}
.right:after {
    content: '';
    position: absolute;
    top: 0; 
    border-top: 15px solid transparent;
    border-left: 15px solid var(--my-color-var);
    border-bottom: 15px solid transparent;
    width: 0;
}

.right_s {    
    position: relative;
    outline: none;
    width: -webkit-calc(100% - 15px);
    width:    -moz-calc(100% - 15px);
    width:         calc(100% - 15px);
    left: -6px;
    z-index: 2;
    transition: 0.3s;
}
.right_s:after {
    content: '';
    position: absolute;
    top: 0; 
    border-top: 3px solid transparent;
    border-left: 6px solid var(--my-color-var);
    border-bottom: 3px solid transparent;
    width: 0;
}


.right_m {    
    position: relative;
    outline: none;
    width: -webkit-calc(100% - 15px);
    width:    -moz-calc(100% - 15px);
    width:         calc(100% - 15px);
    left: -8px;
    z-index: 2;
    transition: 0.3s;
}
.right_m:after {
    content: '';
    position: absolute;
    top: 0; 
    border-top: 8px solid transparent;
    border-left: 8px solid var(--my-color-var);
    border-bottom: 8px solid transparent;
    width: 0;
}

.hidden {
    display: none;
}

.info_box {
    position: fixed;
    top: 0;
    right: -350px;
    width: 300px;
    height: 100%;
    text-align: center;
    z-index: 100;
    background-color: rgba(194, 192, 192, 0.9);
}

.info_box p {
    margin-right: auto;
    margin-left: auto;
    text-align: left;
    width: 80%;
    overflow: visible;
    word-break: break-all;   
}

.info_box textarea {
    margin-right: auto;
    margin-left: auto;
    text-align: left;
    width: 80%;
    height: 300px;
     overflow: visible;
     resize: vertical;
}

.btn {
    position: relative;
    width: 100px;
    height: 30px;
    background-color: rgb(131, 131, 131);
    display: inline-block;
    border: none;
    outline: none;
    margin-top: 15px;
    transition: 0.3s;
    color: white;
    cursor: pointer;  
}

.btn:hover {
    background-color: rgb(70, 181, 201);
    transition: 0.3s;
}

.btn:active {
    background-color: rgb(131, 131, 131);
    transition: 0.3s;
}

.exit {
    position: absolute;
    height: 25px;
    width: 25px;
    outline: none;
    top: 10px;
    right: 10px;
    border: none;
    border-radius: 100%;
    background-color: rgb(131, 131, 131);
    text-align: center;
    transition: 0.3s;
    color: white;
    cursor: pointer;  
}

.exit:hover {
    background-color: red;
    transition: 0.3s;
}

.exit:active {
    background-color: rgb(134, 0, 0);
    transition: 0.3s;
}

#statusmessage{
    position:fixed;
    top:0;
    padding:0 5px;
    line-height:25px;
    background-color:#68d1f1;
    margin-top:-25px;
    font-weight:bold;
    right:0;
    text-align: center;
    z-index: 500;
    width: 300px;

}

</style>

<script>
    $(document).on("mouseenter", ".all", function() {
        var id = $(this).attr("id");
        var arrow = "#" + id + "_arrow"

        $(this).css({ opacity: 1 });
        $(this).css({ 'z-index': 10 });        

        $(arrow).css({ 'z-index': 10 }); 
        $(arrow).css({ opacity: 1 });
    });

    $(document).on("mouseenter", ".line", function() {
        $(this).css({ opacity: 0.8 });
        $(this).css({ 'z-index': 10 });
    });    

    $(document).on("mouseleave", ".all", function() {
        var id = $(this).attr("id");
        var arrow = "#" + id + "_arrow"

        var opc = $("#opacity_slider").val();
        val = opc / 100

        $(this).css({ opacity: val });
        $(this).css({ 'z-index': 2 });

        $(arrow).css({ 'z-index': 2 }); 
        $(arrow).css({ opacity: val });
    });

    $(document).on("mouseleave", ".line", function() {
        var opc = $("#line_opacity_slider").val();
        val = opc / 100

        $(this).css({ opacity: val });
         $(this).css({ 'z-index': 1 });
    });

    function reveal() {
        $(document).on("click", ".all", function() {            

            var id = $(this).attr("id");
            var info = "#" + id + "_info"
            if ($(info).hasClass("hidden")) {
                $(".info_box").addClass("hidden");
                $(info).removeClass("hidden");
                $(info).animate({ right: '0'}, 'slow');  
            } else {
                $(".info_box").animate({ right: '-350'}, 'slow').addClass("hidden");
            }         
        });  
    }
    reveal();

function copy_to_clipboard(id) {
    document.getElementById(id).select();
    document.execCommand('copy');
}

$(document).on("click", ".copy", function() {
    $('#statusmessage').text('Sequence copied to clipboard!').animate({'margin-top':0},200);
    setTimeout( function(){
        $('#statusmessage').animate({'margin-top':-25},200);
    }, 3*1000);
});

$(document).on("click", ".exit", function() {
    $(this).closest(".info_box").animate({ right: '-350'}, 'slow').addClass("hidden");
});

$(document).on("change", "#opacity_slider", function() {
    var opc = $(this).val();
    val = opc / 100
    $(".left, .right, .left_m, .right_m, .left_s, .right_s, .all").css({ opacity: val });
});

$(document).on("change", "#line_opacity_slider", function() {
    var opc = $(this).val();
    val = opc / 100
    $(".line").css({ opacity: val });
});

$(document).on("click", "#small", function() {
    $(".size").removeClass("active");
    $(this).addClass("active");

    $(".R").css({ marginRight : "6px" })
    $(".L").css({ marginLeft : "6px" })

    $(".all").css({ height : "6px" })
    $(".line").css({ height : "2px" })

    $(".AR").removeClass("right");
    $(".AR").removeClass("right_m");
    $(".AR").removeClass("right_s");
    $(".AR").addClass("right_s")

    $(".AL").removeClass("left");
    $(".AL").removeClass("left_m");
    $(".AL").removeClass("left_s");
    $(".AL").addClass("left_s")

    $(".tree").addClass('hidden');
    $(".s_tree").removeClass('hidden')

});

$(document).on("click", "#medium", function() {
    $(".size").removeClass("active");
    $(this).addClass("active");

    $(".R").css({ marginRight : "8px" })
    $(".L").css({ marginLeft : "8px" })

    $(".all").css({ height : "16px" })
    $(".line").css({ height : "8px" })

    $(".AL").removeClass("left");
    $(".AL").removeClass("left_m");
    $(".AL").removeClass("left_s");
    $(".AL").addClass("left_m")

    $(".AR").removeClass("right");
    $(".AR").removeClass("right_m");
    $(".AR").removeClass("right_s");
    $(".AR").addClass("right_m")

    $(".tree").addClass('hidden');
    $(".m_tree").removeClass('hidden')

});

$(document).on("click", "#large", function() {
    $(".size").removeClass("active");
    $(this).addClass("active");

    $(".R").css({ marginRight : "15px" })
    $(".L").css({ marginLeft : "15px" })

    $(".all").css({ height : "30px" })
    $(".line").css({ height : "10px" })

    $(".AL").removeClass("left");
    $(".AL").removeClass("left_m");
    $(".AL").removeClass("left_s");
    $(".AL").addClass("left")

    $(".AR").removeClass("right");
    $(".AR").removeClass("right_m");
    $(".AR").removeClass("right_s");
    $(".AR").addClass("right");

    $(".tree").addClass('hidden');
    $(".l_tree").removeClass('hidden')

});


$(document).on("change", "#scale", function() {
    var scale = $(this).val();
    if (scale == 0) {
        $(".sequences").css({width : "auto"});
    } else {
        $(".sequences").css({width : scale+"px"});
    }
    
});

$(document).ready(function() {
    var w = $( window ).width();
    var p = w - 300;
    $(".sequences").css({width : p});
    $("#scale").val(p);
});

$(document).on("click", "#scale_reset", function() {
    var w = $( window ).width();
    var p = w - 300;
    $(".sequences").css({width : p});  
    $("#scale").val(p);  
});


</script>

<body>

<div class="navigation">
    <button id="small" class="size">Small</button>
    <button id="medium" class="size">Medium</button>
    <button id="large" class="size active">Large</button>
    <input id="opacity_slider" type="range" min="0" max="100" value="70">
    <input id="line_opacity_slider" type="range" min="0" max="100" value="20">
    <input id="scale" type="range" min="0" max="10000">
    <button id="scale_reset">Reset Scale</button>
</div>

<div class="space"></div>

<div id="statusmessage"></div>
<div class="large_grid_container">
    <div class="phylo" style="grid-column-start: 1; grid-column-end: 2;grid-row-start: 1;grid-row-end: 2;">
        <div class="phylo_div">"""
        
html2 = """
        </div>  
    </div>
"""

html3 = """
</div>
</body>
</html>"""

output = html + tree_string + html2 + string + html3

if os.path.exists("test.html"):
    os.remove("test.html")
f = open("test.html", "a")
f.write(output)
f.write('\n')
f.close()
