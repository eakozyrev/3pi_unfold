from array import array

file_read = "data/fit.f"

def parsing(text_start):
    data_points = array('d')
    file = open(file_read)
    str = file.readlines()
    i = 0
    for line in str:
        line = (line.lstrip()).rstrip()
        if line==text_start:
            while str[i+1].find('&') >= 0:
                line1 = str[i+1]
                line1 = (line1.lstrip()).rstrip()
                line1 = (line1.strip('&')).lstrip().strip('/')
                mass = line1.split(',')
                for el in mass:
                    if el!='': data_points.append(float(el))
                i = i+1
            break
        i = i + 1
#    for num in data_points:
	#print(num)
    #print(len(data_points))
    return data_points




#parsing("data ev/")
#print("stop")
#parsing("data edge/")


def save_ev():
    data = parsing("data ev/") + parsing("data evh/")
    print(data,len(data))
    with open("data/ev.dat","w") as file:
        for el in data:
            file.write(str(el))
            file.write('\n')



save_ev()
