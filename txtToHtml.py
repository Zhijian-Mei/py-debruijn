
if __name__ == '__main__':
    froot = 'avastin_5-10mer_0.6_2'
    path = f'{froot}/{froot}_TemplateMatchReport_pre.txt'
    infile = open(path,'r+')
    lines = infile.readlines()
    for line in lines:
        print(len(line.rstrip()))
    quit()
    infile.close()
    message = '''<!DOCTYPE html>
    <body>
    
    '''
    with open(f'{froot}/{froot}_TemplateMatchReport.html','w+') as f:
        quit()
    quit()
    for line in lines:
        print(line.rstrip())