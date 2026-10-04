python3 -c "
from bs4 import BeautifulSoup

with open('cata6.9_2026.html', 'r', encoding='utf-8') as f:
    soup = BeautifulSoup(f, 'html.parser')

rows = soup.find_all('tr', role='row')
lines = []
for row in rows:
    time_td = row.find('td', class_='time')
    if not time_td:
        continue
    a_tag = time_td.find('a')
    if not a_tag:
        continue
    href = a_tag.get('href', '')
    event_id = href.rstrip('/').split('/')[-1]
    
    time_str = a_tag.get_text(strip=True)
    parts = time_str.split()
    if len(parts) != 2:
        continue
    date_part, clock_part = parts[0], parts[1]
    tt = f'{date_part}-{clock_part.replace(\":\", \"\")}'
    
    lat_td = row.find('td', class_='latitude')
    lon_td = row.find('td', class_='longitude')
    dep_td = row.find('td', class_='depth')
    mag_td = row.find('td', class_='magnitude')
    
    if not (lat_td and lon_td and dep_td and mag_td):
        continue
        
    lat = lat_td.get_text(strip=True)
    lon = lon_td.get_text(strip=True)
    dep = dep_td.get_text(strip=True)
    mag = mag_td.get_text(strip=True)
    
    iso_time = f'{date_part}T{clock_part}'
    lines.append(f'{tt} {mag} {event_id} {lon} {lat} {dep} {iso_time}')

with open('weblist', 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines) + ('\n' if lines else ''))
"

cat weblist | gawk '{print "mkdir -p "$1"_"$2"_"$3}' | sh
cat weblist | gawk '{print "echo "$4,$5,$6,$7" > "$1"_"$2"_"$3"/loladep.dat"}' | sh

cat weblist | gawk '{print "sh sub_fetchdata.sh "$1,$2,$3}' | sh
cat weblist | gawk '{print "sh sub_processdata.sh "$1,$2,$3}' | sh

# need one more step to merge data to ../IRISloc/
