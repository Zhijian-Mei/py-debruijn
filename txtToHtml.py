
if __name__ == '__main__':
    froot = 'avastin_5-10mer_0.6_2'
    path = f'{froot}/{froot}_TemplateMatchReport.txt'
    infile = open(path,'r+')
    lines = infile.readlines()
    infile.close()
    message = '''<!DOCTYPE html>
    <body>
    
    </body>
    '''
    with open(f'{froot}/{froot}_TemplateMatchReport.html','w+') as f:
        f.write(message)
    quit()
    for line in lines:
        print(line.rstrip())