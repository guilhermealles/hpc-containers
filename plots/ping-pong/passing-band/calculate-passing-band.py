import sys

file = open(sys.argv[1], 'r')
lines = file.read().split('\n')

passing_band_text = ''

for line in lines:
    if line == '':
        break
    size_bytes = int(line.split(' ')[0])
    time = float(line.split(' ')[1])/1000

    kilobytes_per_second = (float(size_bytes)/(1024))/time
    current_line = str(size_bytes) + ' ' + str(kilobytes_per_second) + '\n'
    passing_band_text = passing_band_text + current_line

print(passing_band_text)

