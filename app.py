import streamlit as st
from datetime import datetime, timedelta
from math import sin, cos, tan, atan2, degrees, radians
from astropy.time import Time
from astropy.coordinates import get_body, solar_system_ephemeris, GeocentricTrueEcliptic
from collections import defaultdict
import pandas as pd
from geopy.geocoders import Nominatim
from geopy.extra.rate_limiter import RateLimiter
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import io
# NEW: timezone detection for "Current City"
from timezonefinder import TimezoneFinder
import pytz
# ---- Matplotlib defaults (crisp + thin) ----
plt.rcParams.update({"figure.dpi": 300, "savefig.dpi": 300, "lines.linewidth": 0.28})
# ---- Constants ----
cities_fallback = {
    'Chennai': {'lat': 13.08, 'lon': 80.27}, 'Mumbai': {'lat': 19.07, 'lon': 72.88},
    'Delhi': {'lat': 28.61, 'lon': 77.23}, 'Bangalore': {'lat': 12.97, 'lon': 77.59},
    'Kolkata': {'lat': 22.57, 'lon': 88.36}, 'Hyderabad': {'lat': 17.39, 'lon': 78.49},
}
sign_names = ['Aries','Taurus','Gemini','Cancer','Leo','Virgo','Libra','Scorpio','Sagittarius','Capricorn','Aquarius','Pisces']
lords_full = ['Ketu','Venus','Sun','Moon','Mars','Rahu','Jupiter','Saturn','Mercury']
lords_short = ['Ke','Ve','Su','Mo','Ma','Ra','Ju','Sa','Me']
nak_names = ['Ashwini','Bharani','Krittika','Rohini','Mrigashira','Ardra','Punarvasu','Pushya','Ashlesha',
             'Magha','Purva Phalguni','Uttara Phalguni','Hasta','Chitra','Swati','Vishakha','Anuradha',
             'Jyeshta','Mula','Purva Ashadha','Uttara Ashadha','Shravana','Dhanishta','Shatabhisha',
             'Purva Bhadrapada','Uttara Bhadrapada','Revati']
years = [7, 20, 6, 10, 7, 18, 16, 19, 17] * 3
sign_lords = ['Mars','Venus','Mercury','Moon','Sun','Mercury','Venus','Mars','Jupiter','Saturn','Saturn','Jupiter']
sthana_bala_dict = {
    'Sun': [100,90,80,70,60,50,40,50,60,70,80,90],
    'Moon': [90,100,90,80,70,60,60,50,70,70,70,90],
    'Jupiter': [60,60,70,100,90,60,75,60,80,40,50,80],
    'Venus': [60,70,60,50,40,35,80,50,60,80,70,100],
    'Mercury': [40,60,70,45,60,100,60,45,55,50,45,35],
    'Mars': [80,70,45,35,60,45,50,60,60,100,90,60],
    'Saturn': [35,50,60,70,80,60,100,90,50,60,80,50],
    'Rahu': [100]*12,
    'Ketu': [100]*12
}
capacity_dict = {
    'Saturn': 100,
    'Mars': 50,
    'Sun': 100,
    'Jupiter': 100,
    'Venus': 50,
    'Mercury': 30,
    'Moon': 100,
    'Rahu': 100,
    'Ketu': 50
}
good_capacity_dict = {
    'Saturn': 0,
    'Mars': 75,
    'Sun': 50,
    'Jupiter': 100,
    'Venus': 100,
    'Mercury': 100,
    'Rahu': 0,
    'Ketu': 100
}
bad_capacity_dict = {
    'Saturn': 100,
    'Mars': 25,
    'Sun': 50,
    'Jupiter': 0,
    'Venus': 0,
    'Mercury': 0,
    'Rahu': 100,
    'Ketu': 0
}
shukla_good = [100, 9, 16, 23, 30, 37, 44, 51, 58, 65, 72, 79, 86, 93, 100]
shukla_bad = [0] * 15
krishna_good = [93, 86, 79, 72, 65, 58, 51, 44, 37, 30, 23, 16, 9, 2, 0]
krishna_bad = [7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 84, 91, 98, 100]
order_dict = {'Jupiter':1, 'Venus':2, 'Mercury':3, 'Sun':4, 'Mars':5, 'Saturn':6, 'Rahu':7, 'Ketu':8}
mix_dict = {0:100,1:100,2:100,3:95,4:90,5:85,6:80,7:75,8:70,9:65,10:60,11:55,12:50,13:45,14:40,15:35,16:30,17:25,18:20,19:15,20:10,21:5,22:0}
# ---- Astro helpers ----
def get_lahiri_ayanamsa(year):
    base = 23.853; rate = 50.2388/3600.0
    return (base + (year - 2000) * rate) % 360
def get_obliquity(d):
    T = d/36525.0
    return ((((-4.34e-8*T - 5.76e-7)*T + 0.0020034)*T - 1.831e-4)*T - 46.836769)*T/3600 + 23.4392794444444
def get_gmst(d):
    T = d/36525.0
    return (67310.54841 + (3155760000 + 8640184.812866)*T + 0.093104*T**2 - 6.2e-6*T**3)/3600 % 24
def get_ascendant(jd, lat, lon):
    d = jd - 2451545.0
    oer = radians(get_obliquity(d))
    lst = (get_gmst(d) + lon/15.0) % 24
    lstr = radians(lst*15.0)
    sin_asc = cos(lstr)
    cos_asc = -(sin(lstr)*cos(oer) + tan(radians(lat))*sin(oer))
    return degrees(atan2(sin_asc, cos_asc)) % 360
def get_sidereal_lon(tlon, ayan): return (tlon - ayan) % 360
def get_sign(lon): return sign_names[int(lon/30)]
def get_house(lon, lagna_lon): return (int(lon/30) - int(lagna_lon/30)) % 12 + 1
def get_nakshatra_details(lon):
    dnak = 360/27
    idx = int(lon // dnak) % 27
    pos = lon % dnak
    pada = int((pos/(dnak/4))) + 1
    star = idx % 9
    sub = (star + int((pos/dnak)*9)) % 9
    return nak_names[idx], pada, lords_short[star], lords_short[sub]
def generate_vimshottari_dasa(moon_lon):
    nak = int(moon_lon * 27 / 360)
    lord_idx = nak % 9
    y = years[lord_idx]
    pos_in_nak = moon_lon % (360/27)
    fraction = pos_in_nak / (360/27)
    return lord_idx, y * (1 - fraction)
def generate_periods(start_date, lord_idx, total_years, level='dasa', max_depth=3):
    periods, remaining, i, current = [], total_years, lord_idx, start_date
    depth_map = {'dasa':0,'bhukti':1,'anthara':2,'sukshma':3,'prana':4,'sub_prana':5}
    next_level = {0:'bhukti',1:'anthara',2:'sukshma',3:'prana',4:'sub_prana',5:None}
    depth = depth_map.get(level,0)
    while remaining > 0:
        lord = lords_full[i]; y_full = years[i]
        y = min((y_full/120)*total_years, remaining)
        end = current + timedelta(days=y*365.25)
        subs = generate_periods(current, i, y, next_level[depth], max_depth) if (depth < max_depth-1 and next_level[depth]) else []
        periods.append((lord, current, end, subs))
        remaining -= y; current = end; i = (i+1) % 9
    return periods
def filter_from_birth(periods, birth_dt):
    out = []
    for lord, start, end, sub in periods:
        if end > birth_dt:
            out.append((lord, max(start, birth_dt), end, filter_from_birth(sub, birth_dt) if sub else []))
    return out
def duration_str(delta, level='dasa'):
    days = delta.total_seconds()/86400
    if days < 1 and level in ['sukshma','prana','sub_prana']:
        hrs = days*24; h = int(hrs); m = int((hrs - h)*60)
        return "Less than 1 minute" if h==0 and m==0 else f"{h}h {m}m"
    y = int(days/365.25); rem = days % 365.25
    m = int(rem/30.4375); d = int(rem % 30.4375)
    return "Less than 1 day" if y+m+d==0 else f"{y}y {m}m {d}d"
def calculate_dig_bala(planet, lon, lagna):
    east = lagna % 360
    north = (lagna + 90) % 360
    west = (lagna + 180) % 360
    south = (lagna + 270) % 360
    if planet.lower() in ['sun', 'mars']:
        D = north
    elif planet.lower() in ['moon', 'venus']:
        D = south
    elif planet.lower() in ['mercury', 'jupiter', 'ketu']:
        D = west
    elif planet.lower() in ['saturn', 'rahu']:
        D = east
    else:
        return None
    diff = abs(lon - D)
    ang_dist = min(diff, 360 - diff)
    virupas = ang_dist / 3
    percentage = (virupas / 60) * 100
    return round(percentage, 2)
def compute_chart(name, date_obj, time_str, lat, lon, tz_offset, max_depth):
    # parse time
    try:
        hour, minute = map(int, time_str.split(':'))
        if not (0<=hour<=23 and 0<=minute<=59): raise ValueError
    except:
        raise ValueError("Time must be in HH:MM format (24-hour)")
    local_dt = datetime.combine(date_obj, datetime.min.time().replace(hour=hour, minute=minute))
    utc_dt = local_dt - timedelta(hours=tz_offset)
    t = Time(utc_dt); jd = t.jd; ayan = get_lahiri_ayanamsa(utc_dt.year)
    with solar_system_ephemeris.set('builtin'):
        lon_trop = {}
        for nm in ['sun','moon','mercury','venus','mars','jupiter','saturn']:
            ecl = get_body(nm, t).transform_to(GeocentricTrueEcliptic()); lon_trop[nm] = ecl.lon.deg
    d = jd - 2451545.0; T = d/36525.0
    omega = (125.04452 - 1934.136261*T + 0.0020708*T**2 + T**3/450000) % 360
    lon_trop['rahu'] = omega; lon_trop['ketu'] = (omega + 180) % 360
    lon_sid = {p: get_sidereal_lon(lon_trop[p], ayan) for p in lon_trop}
    lagna_sid = get_sidereal_lon(get_ascendant(jd, lat, lon), ayan)
    sun_lon = lon_sid['sun']
    moon_lon = lon_sid['moon']
    diff = (moon_lon - sun_lon) % 360
    tithi_fraction = diff / 12
    tithi = int(tithi_fraction) + 1
    if tithi > 30:
        tithi = 30
    if diff < 180:
        paksha = 'Krishna'  # Moon before Sun, waning
    else:
        paksha = 'Shukla'  # Moon after Sun, waxing
    if paksha == 'Shukla':
        tithi_idx = tithi - 1 if tithi <=15 else tithi -16
    else:
        tithi_idx = tithi -1 if tithi <=15 else tithi -16
    # rasi houses first for conjunctions
    house_planets_rasi = defaultdict(list)
    positions = {**lon_sid, 'asc': lagna_sid}
    for p, L in positions.items():
        house_planets_rasi[get_house(L, lagna_sid)].append(p.capitalize() if p != 'asc' else 'Asc')
    # planets table
    rows = []
    planet_data = {}
    consumed_notes = {}
    asc_deg = lagna_sid % 360; asc_sign = get_sign(asc_deg)
    a_nak, a_pada, a_ld, a_sl = get_nakshatra_details(asc_deg)
    dig_bala_asc = calculate_dig_bala('asc', asc_deg, lagna_sid)
    rows.append(['Asc', f"{asc_deg:.2f}", asc_sign, a_nak, a_pada, f"{a_ld}/{a_sl}", f"{dig_bala_asc}%" if dig_bala_asc is not None else '', '', '', '', '', ''])
    for p in ['sun','moon','mars','mercury','jupiter','venus','saturn','rahu','ketu']:
        L = lon_sid[p]; sign = get_sign(L); nak, pada, ld, sl = get_nakshatra_details(L)
        dig_bala = calculate_dig_bala(p, L, lagna_sid)
        planet_cap = p.capitalize()
        sthana = sthana_bala_dict.get(planet_cap, [0]*12)[sign_names.index(sign)]
        capacity = capacity_dict.get(planet_cap, None)
        volume = (capacity * sthana / 100.0) if capacity is not None else ''
        cap_good = volume
        cap_bad = volume
        if planet_cap == 'Rahu':
            good_volume = 0.0
            bad_volume = volume
        elif planet_cap == 'Ketu':
            good_volume = 50.0
            bad_volume = 0.0
        else:
            if planet_cap == 'Moon':
                if paksha == 'Shukla':
                    good_capacity = shukla_good[tithi_idx]
                    bad_capacity = shukla_bad[tithi_idx]
                else:
                    good_capacity = krishna_good[tithi_idx]
                    bad_capacity = krishna_bad[tithi_idx]
            else:
                good_capacity = good_capacity_dict.get(planet_cap, None)
                bad_capacity = bad_capacity_dict.get(planet_cap, None)
            good_volume = (volume * (good_capacity / 100.0)) if good_capacity is not None and volume != '' else ''
            bad_volume = (volume * (bad_capacity / 100.0)) if bad_capacity is not None and volume != '' else ''
        planet_data[planet_cap] = {'sthana': sthana, 'volume': volume, 'good_volume': good_volume, 'bad_volume': bad_volume, 'cap_good': cap_good, 'cap_bad': cap_bad, 'dig_bala': dig_bala, 'L': L, 'sign': sign, 'nak': nak, 'pada': pada, 'ld_sl': f"{ld}/{sl}"}
        consumed_notes[planet_cap] = {'good': [], 'bad': []}
    # Adjust for conjunctions
    for h in range(1, 13):
        house_planets = [p for p in house_planets_rasi[h] if p != 'Asc']
        if len(house_planets) > 1:
            # Special for Ketu
            if 'Ketu' in house_planets and ('Sun' in house_planets or 'Moon' in house_planets):
                planet_data['Ketu']['bad_volume'] = planet_data['Ketu']['volume']
                planet_data['Ketu']['good_volume'] = 0.0
            # General grab for any with bad_volume >0
            general_grabbers = [p for p in house_planets if planet_data[p]['bad_volume'] > 0]
            general_grabbers.sort(key=lambda p: -order_dict.get(p, 0))  # higher order first
            for grabber in general_grabbers:
                for grabbed in [p for p in house_planets if p != grabber and planet_data[p]['good_volume'] > 0]:
                    if grabber == 'Ketu' and grabbed not in ['Sun', 'Moon']:
                        continue
                    if grabber == 'Jupiter' and grabbed != 'Moon':
                        continue
                    deg_diff = int(abs(planet_data[grabber]['L'] - planet_data[grabbed]['L']))
                    mix = mix_dict.get(min(deg_diff, 22), 0) / 100.0
                    available_grab = planet_data[grabbed]['good_volume'] * mix
                    room = planet_data[grabber]['cap_good'] - planet_data[grabber]['good_volume']
                    grab_amount = min(available_grab, room)
                    if grab_amount > 0:
                        planet_data[grabbed]['good_volume'] -= grab_amount
                        planet_data[grabber]['good_volume'] += grab_amount
                        consumed_notes[grabber]['good'].append(f"{grab_amount:.2f} from {grabbed}")
                        consumed_notes[grabbed]['good'].append(f"{-grab_amount:.2f} to {grabber}")
            # For Ketu exchange bad in full
            if 'Ketu' in house_planets:
                ketu_bad = planet_data['Ketu']['bad_volume']
                grabbed_list = [p for p in house_planets if p in ['Sun', 'Moon']]
                if grabbed_list:
                    add_per = ketu_bad / len(grabbed_list)
                    for grabbed in grabbed_list:
                        room_bad = planet_data[grabbed]['cap_bad'] - planet_data[grabbed]['bad_volume']
                        add = min(add_per, room_bad)
                        planet_data[grabbed]['bad_volume'] += add
                        consumed_notes[grabbed]['bad'].append(f"{add:.2f} from Ketu")
                        consumed_notes['Ketu']['bad'].append(f"{-add:.2f} to {grabbed}")
            # For good planets exchange
            remaining_bad = sum(planet_data[p]['bad_volume'] for p in house_planets)
            if remaining_bad == 0:
                good_planets = [p for p in house_planets if planet_data[p]['good_volume'] > 0]
                if len(good_planets) > 1:
                    original_goods = {p: planet_data[p]['good_volume'] for p in good_planets}
                    total_good = sum(original_goods.values())
                    avg_good = total_good / len(good_planets)
                    losers = [p for p in good_planets if original_goods[p] > avg_good]
                    total_loss = sum(original_goods[p] - avg_good for p in losers)
                    for p in good_planets:
                        planet_data[p]['good_volume'] = min(avg_good, planet_data[p]['cap_good'] - planet_data[p]['bad_volume'])
                    for p in good_planets:
                        diff = planet_data[p]['good_volume'] - original_goods[p]
                        if diff > 0:
                            for loser in losers:
                                proportion = (original_goods[loser] - avg_good) / total_loss if total_loss > 0 else 0
                                add_from = diff * proportion
                                consumed_notes[p]['good'].append(f"{add_from:.2f} from {loser}")
                        elif diff < 0:
                            consumed_notes[p]['good'].append(f"{diff:.2f} to group")
            # For good planets that lost good, share back based on degree gap
            lost_good_planets = [p for p in house_planets if planet_data[p]['bad_volume'] == 0 and planet_data[p]['good_volume'] < planet_data[p]['volume'] * (good_capacity_dict.get(p, 100) / 100)]
            for lost in lost_good_planets:
                share_from = [p for p in house_planets if p != lost and planet_data[p]['good_volume'] > 0]
                for sharer in share_from:
                    deg_diff = int(abs(planet_data[lost]['L'] - planet_data[sharer]['L']))
                    mix = mix_dict.get(min(deg_diff, 22), 0) / 100.0
                    available_share = planet_data[sharer]['good_volume'] * mix
                    needed = planet_data[lost]['volume'] * (good_capacity_dict.get(lost, 100) / 100) - planet_data[lost]['good_volume']
                    share_amount = min(available_share, needed)
                    if share_amount > 0:
                        planet_data[sharer]['good_volume'] -= share_amount
                        planet_data[lost]['good_volume'] += share_amount
                        consumed_notes[lost]['good'].append(f"{share_amount:.2f} shared from {sharer}")
                        consumed_notes[sharer]['good'].append(f"{-share_amount:.2f} shared to {lost}")
    # Build rows with adjusted
    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        data = planet_data[p]
        notes = []
        if consumed_notes[p]['good']:
            notes.append(f"Good: {'; '.join(consumed_notes[p]['good'])}")
        if consumed_notes[p]['bad']:
            notes.append(f"Bad: {'; '.join(consumed_notes[p]['bad'])}")
        note_str = '; '.join(notes)
        rows.append([p, f"{data['L']:.2f}", data['sign'], data['nak'], data['pada'], data['ld_sl'], f"{data['dig_bala']}%" if data['dig_bala'] is not None else '', f"{data['sthana']}%", f"{data['volume']:.2f}" if isinstance(data['volume'], float) else '', f"{data['good_volume']:.2f}" if isinstance(data['good_volume'], float) else '', f"{data['bad_volume']:.2f}" if isinstance(data['bad_volume'], float) else '', note_str])
    df_planets = pd.DataFrame(rows, columns=['Planet','Deg','Sign','Nakshatra','Pada','Ld/SL','Dig Bala (%)','Sthana Bala (%)','Volume','Good Volume','Bad Volume','Consumed Notes'])
    # df_rasi
    df_rasi = pd.DataFrame([[f"House {h}", get_sign((lagna_sid+(h-1)*30)%360),
                             ', '.join(sorted(house_planets_rasi[h])) if house_planets_rasi[h] else 'Empty']
                            for h in range(1,13)], columns=['House','Sign','Planets'])
    # navamsa
    nav_lagna = (lagna_sid*9) % 360
    house_planets_nav = defaultdict(list)
    for p,L in lon_sid.items():
        nav_lon = (L*9) % 360
        nav_h = (int(nav_lon/30) - int(nav_lagna/30)) % 12 + 1
        house_planets_nav[nav_h].append(p.capitalize())
    df_nav = pd.DataFrame([[f"House {h}", get_sign((nav_lagna+(h-1)*30)%360),
                             ', '.join(sorted(house_planets_nav[h])) if house_planets_nav[h] else 'Empty']
                            for h in range(1,13)], columns=['House','Sign','Planets'])
    # aspects table
    lagna_sign = get_sign(lagna_sid)
    aspects_dict = {'Sun':[7],'Moon':[7],'Mars':[4,7,8],'Mercury':[7],'Jupiter':[5,7,9],'Venus':[7],'Saturn':[3,7,10]}
    planet_to_house = {p.capitalize(): get_house(lon_sid[p], lagna_sid) for p in lon_sid}
    house_status = []
    for h in range(1,13):
        lord = sign_lords[(sign_names.index(lagna_sign)+(h-1))%12]
        lord_house = planet_to_house[lord]
        asp = []
        for planet, offs in aspects_dict.items():
            if planet in planet_to_house:
                ph = planet_to_house[planet]
                for off in offs:
                    if ((ph-1+(off-1))%12)+1 == h: asp.append(planet)
        house_status.append([f"House {h}",
                             ', '.join(sorted(house_planets_rasi[h])) if house_planets_rasi[h] else 'Empty',
                             ', '.join(asp) if asp else 'None', lord, f"House {lord_house}"])
    df_house_status = pd.DataFrame(house_status, columns=['House','Planets','Aspects from','Lord','Lord in'])
    # dasa tree
    moon_lon = lon_sid['moon']
    idx, bal = generate_vimshottari_dasa(moon_lon)
    full_first = years[idx]; passed = full_first - bal
    dasa_start = utc_dt - timedelta(days=passed*365.25)
    dasa = generate_periods(dasa_start, idx, 120, 'dasa', max_depth)
    dasa_filtered = filter_from_birth(dasa, utc_dt)
    depth_map = {1:'Dasa only',2:'Dasa + Bhukti',3:'Dasa + Bhukti + Anthara',
                 4:'Dasa + Bhukti + Anthara + Sukshma',5:'Dasa + Bhukti + Anthara + Sukshma + Prana',
                 6:'Dasa + Bhukti + Anthara + Sukshma + Prana + Sub-Prana'}
    return {
        'name': name, 'df_planets': df_planets, 'df_rasi': df_rasi, 'df_nav': df_nav,
        'df_house_status': df_house_status, 'dasa_periods_filtered': dasa_filtered,
        'lagna_sid': lagna_sid, 'nav_lagna': nav_lagna, 'lagna_sign': lagna_sign,
        'nav_lagna_sign': get_sign(nav_lagna), 'moon_rasi': get_sign(moon_lon),
        'moon_nakshatra': get_nakshatra_details(moon_lon)[0], 'moon_pada': get_nakshatra_details(moon_lon)[1],
        'selected_depth': depth_map[max_depth], 'utc_dt': utc_dt, 'max_depth': max_depth,
        'house_to_planets_rasi': house_planets_rasi, 'house_to_planets_nav': house_planets_nav
    }
# ---- South Indian plotter ----
def plot_south_indian_style(ax, house_to_planets, lagna_sign, title):
    sign_positions = {'Pisces':(0,3),'Aries':(1,3),'Taurus':(2,3),'Gemini':(3,3),
                      'Cancer':(3,2),'Leo':(3,1),'Virgo':(3,0),
                      'Libra':(2,0),'Scorpio':(1,0),'Sagittarius':(0,0),
                      'Capricorn':(0,1),'Aquarius':(0,2)}
    lagna_idx = sign_names.index(lagna_sign)
    house_for_sign = {s: ((i - lagna_idx) % 12) + 1 for i, s in enumerate(sign_names)}
    box_w, box_h, spacing, pad = 0.46, 0.46, 0.52, 0.02
    top_pad_extra = 0.020
    line_h_min, line_h_max = 0.042, 0.058
    planet_font = 2.45
    for sign,(gx,gy) in sign_positions.items():
        h = house_for_sign[sign]
        planets = sorted(house_to_planets.get(h,[]))
        x = gx*spacing + 0.22; y = (3-gy)*spacing + 0.22
        ax.add_patch(FancyBboxPatch((x,y), box_w, box_h,
                                    boxstyle="round,pad=0.004",
                                    ec="black", fc="#F5F5F5",
                                    alpha=0.92, linewidth=0.32))
        ax.text(x+pad, y+pad, sign[:3], ha='left', va='top', fontsize=2.7)
        if planets:
            avail = box_h - (pad + 0.064)
            n = len(planets)
            line_h = min(line_h_max, max(line_h_min, avail / max(1, n)))
            start_y = y + pad + 0.040 + top_pad_extra
            for i, p in enumerate(planets):
                py = start_y + i*line_h
                ax.text(x + box_w/2, py, p, ha='center', va='top', fontsize=planet_font)
    ax.set_xlim(0,3); ax.set_ylim(0,3); ax.set_aspect('equal'); ax.invert_yaxis()
    ax.set_title(title, fontsize=3.6, fontweight='normal')
    ax.axis('off')
# ---- Streamlit UI ----
st.set_page_config(page_title="Sivapathy Horoscope", layout="wide")
st.markdown("""
<style>
    .stApp { background-color: white; color: #125336; }
    .stTextInput > div > div > input,
    .stSelectbox > div > div > select,
    .stNumberInput > div > div > input {
        background-color: white; color: #125336; border: 1px solid #125336;
    }
    .stButton > button { background-color: #125336; color: white; border: none; padding: 0.5rem 2rem; font-size: 1.1rem; font-weight: 600; }
    .stButton > button:hover { background-color: #0a3d22; box-shadow: 0 4px 6px rgba(0,0,0,0.1); }
    h1, h2, h3 { color: #125336 !important; }
    .summary-box { background-color: #f0f7f4; padding: 1.2rem; border-radius: 10px; border: 2px solid #125336; margin: 1rem 0; }
    .summary-item { font-size: 1.05rem; margin: 0.35rem 0; color: #125336; }
</style>
""", unsafe_allow_html=True)
st.title("Sivapathy Astrology Data Generator")
if 'chart_data' not in st.session_state: st.session_state.chart_data = None
if 'search_results' not in st.session_state: st.session_state.search_results = []
@st.cache_resource
def get_geolocator():
    geolocator = Nominatim(user_agent="vedic_astro_app")
    return RateLimiter(geolocator.geocode, min_delay_seconds=1)
geocode = get_geolocator()
# === Timezone tools for Current City ===
_tf = TimezoneFinder()
def tz_for_latlon(lat: float, lon: float):
    tzname = _tf.timezone_at(lng=lon, lat=lat)
    if not tzname:
        return pytz.UTC
    return pytz.timezone(tzname)
# === Tree walkers for micro-depths ===
_DEPTH_NAME_TO_INT = {'Dasa':1, 'Bhukti':2, 'Anthara':3, 'Sukshma':4, 'Prana':5, 'Sub-Prana':6}
def find_active_path_to_depth(periods, when_utc, target_depth, cur_depth=1):
    for lord, start, end, subs in periods:
        if start <= when_utc < end:
            if cur_depth == target_depth or not subs:
                return [(lord, start, end)]
            sub_path = find_active_path_to_depth(subs, when_utc, target_depth, cur_depth+1)
            return [(lord, start, end)] + (sub_path or [])
    return None
def collect_periods_at_depth(periods, target_depth, cur_depth=1, acc=None):
    if acc is None: acc = []
    for lord, start, end, subs in periods:
        if cur_depth == target_depth or not subs:
            acc.append((lord, start, end))
        else:
            collect_periods_at_depth(subs, target_depth, cur_depth+1, acc)
    return acc
# =========================
# Birth Details (Birth City separate)
# =========================
st.subheader("Birth Details")
name = st.text_input("Name", placeholder="Enter full name")
c1, c2, c3 = st.columns(3)
with c1:
    birth_date = st.date_input("Birth Date", value=datetime.now().date(),
                               min_value=datetime(1900,1,1).date(), max_value=datetime.now().date())
with c2:
    birth_time = st.text_input("Birth Time (HH:MM in 24-hour format)", placeholder="14:30")
with c3:
    tz_offset = st.number_input("Timezone offset at birth (hrs)", value=5.5, step=0.5,
                                help="Offset from UTC at birth (e.g., IST = 5.5)")
use_custom_coords = st.checkbox("Custom birth latitude and longitude?")
if use_custom_coords:
    clat, clon = st.columns(2)
    with clat: lat = st.number_input("Birth Latitude", value=13.08, format="%.4f")
    with clon: lon = st.number_input("Birth Longitude", value=80.27, format="%.4f")
else:
    birth_city_query = st.text_input("Birth City (for birth coordinates)",
                                     placeholder="Start typing birth city name...", key="birth_city_input")
    if birth_city_query and len(birth_city_query) >= 2:
        try:
            locations = geocode(birth_city_query, exactly_one=False, limit=5)
            st.session_state.search_results = [{'display': loc.address, 'lat': loc.latitude, 'lon': loc.longitude, 'address': loc.address} for loc in (locations or [])]
        except:
            st.session_state.search_results = []
    else:
        st.session_state.search_results = []
    if st.session_state.search_results:
        opts = [r['display'] for r in st.session_state.search_results]
        sel = st.selectbox("Select birth location", options=opts, key="birth_location_selector")
        i = opts.index(sel)
        lat = st.session_state.search_results[i]['lat']; lon = st.session_state.search_results[i]['lon']
        st.success(f"Selected birth place: {st.session_state.search_results[i]['address']} (Lat: {lat:.2f}, Lon: {lon:.2f})")
    else:
        city_key = (birth_city_query or "").title()
        if city_key in cities_fallback:
            lat = cities_fallback[city_key]['lat']; lon = cities_fallback[city_key]['lon']
        else:
            lat, lon = 13.08, 80.27
        st.info("Using default birth location: Chennai, India")
max_depth_options = {
    1:'Dasa only',2:'Dasa + Bhukti',3:'Dasa + Bhukti + Anthara',
    4:'Dasa + Bhukti + Anthara + Sukshma',5:'Dasa + Bhukti + Anthara + Sukshma + Prana',
    6:'Dasa + Bhukti + Anthara + Sukshma + Prana + Sub-Prana'
}
selected_depth_str = st.selectbox("Generate up to (depth)", list(max_depth_options.values()), index=2)
max_depth = [k for k,v in max_depth_options.items() if v == selected_depth_str][0]
if st.button("Generate Chart", use_container_width=True):
    if not name:
        st.error("Please enter a name.")
    elif not birth_time:
        st.error("Please enter birth time in HH:MM format.")
    else:
        try:
            with st.spinner("Calculating chart..."):
                st.session_state.chart_data = compute_chart(name, birth_date, birth_time, lat, lon, tz_offset, max_depth)
            st.success("Chart generated successfully!")
            st.rerun()
        except ValueError as e:
            st.error(f"Invalid input: {e}")
        except Exception as e:
            st.error(f"Error generating chart: {e}")
# ---- render helpers
def show_png(fig):
    fig.tight_layout(pad=0.10)
    st.pyplot(fig, use_container_width=False, dpi=300)
# =========================
# Outputs
# =========================
if st.session_state.chart_data:
    cd = st.session_state.chart_data
    st.markdown("---")
    st.markdown(f"""
    <div class="summary-box">
        <h3>Chart Summary</h3>
        <div class="summary-item"><strong>Name:</strong> {cd['name']}</div>
        <div class="summary-item"><strong>Lagna:</strong> {cd['lagna_sign']} ({cd['lagna_sid']:.2f}°)</div>
        <div class="summary-item"><strong>Rasi (Moon Sign):</strong> {cd['moon_rasi']}</div>
        <div class="summary-item"><strong>Nakshatra:</strong> {cd['moon_nakshatra']} (Pada {cd['moon_pada']})</div>
    </div>
    """, unsafe_allow_html=True)
    st.subheader("Planetary Positions")
    st.dataframe(cd['df_planets'], hide_index=True, use_container_width=True)
    # Rasi & Navamsa (South Indian)
    st.subheader("Rasi (D1) & Navamsa (D9) — South Indian")
    col1, col2 = st.columns(2, gap="small")
    size = (1.8, 1.8)
    fig1, ax1 = plt.subplots(figsize=size)
    plot_south_indian_style(ax1, cd['house_to_planets_rasi'], cd['lagna_sign'], 'Rasi Chart (South Indian)')
    show_png(fig1)
    fig2, ax2 = plt.subplots(figsize=size)
    plot_south_indian_style(ax2, cd['house_to_planets_nav'], cd['nav_lagna_sign'], 'Navamsa Chart (South Indian)')
    show_png(fig2)
    st.subheader("House Analysis")
    st.dataframe(cd['df_house_status'], hide_index=True, use_container_width=True)
    st.subheader(f"Vimshottari Dasa ({cd['selected_depth']})")
    dasa_rows = [{'Planet': lord, 'Start': s.strftime('%Y-%m-%d'),
                  'End': e.strftime('%Y-%m-%d'), 'Duration': duration_str(e-s,'dasa')}
                 for lord, s, e, _ in cd['dasa_periods_filtered']]
    st.dataframe(pd.DataFrame(dasa_rows), hide_index=True, use_container_width=True)
    max_depth_sel = cd['max_depth']
    dp = cd['dasa_periods_filtered']
    if max_depth_sel >= 2:
        with st.expander("View Bhuktis (Sub-periods)", expanded=False):
            if dp:
                d_opt = [f"{p[0]} ({p[1].strftime('%Y-%m-%d')} - {p[2].strftime('%Y-%m-%d')})" for p in dp]
                sel = st.selectbox("Select Dasa:", d_opt, key="dasa_select")
                bhuktis = dp[d_opt.index(sel)][3]
                st.dataframe(pd.DataFrame(
                    [{'Planet': l, 'Start': s.strftime('%Y-%m-%d'), 'End': e.strftime('%Y-%m-%d'),
                      'Duration': duration_str(e-s,'bhukti')} for l,s,e,_ in bhuktis]),
                    hide_index=True, use_container_width=True)
                if max_depth_sel >= 3:
                    with st.expander("View Antharas", expanded=False):
                        if bhuktis:
                            b_opt = [f"{p[0]} ({p[1].strftime('%Y-%m-%d')} - {p[2].strftime('%Y-%m-%d')})" for p in bhuktis]
                            selb = st.selectbox("Select Bhukti:", b_opt, key="bhukti_select")
                            antharas = bhuktis[b_opt.index(selb)][3]
                            st.dataframe(pd.DataFrame(
                                [{'Planet': l, 'Start': s.strftime('%Y-%m-%d %H:%M'),
                                  'End': e.strftime('%Y-%m-%d %H:%M'),
                                  'Duration': duration_str(e-s,'anthara')} for l,s,e,_ in antharas]),
                                hide_index=True, use_container_width=True)
                            if max_depth_sel >= 4:
                                with st.expander("View Sukshmas", expanded=False):
                                    if antharas:
                                        a_opt = [f"{p[0]} ({p[1].strftime('%Y-%m-%d %H:%M')} - {p[2].strftime('%Y-%m-%d %H:%M')})" for p in antharas]
                                        sela = st.selectbox("Select Anthara:", a_opt, key="anthara_select")
                                        sukshmas = antharas[a_opt.index(sela)][3]
                                        st.dataframe(pd.DataFrame(
                                            [{'Planet': l, 'Start': s.strftime('%Y-%m-%d %H:%M'),
                                              'End': e.strftime('%Y-%m-%d %H:%M'),
                                              'Duration': duration_str(e-s,'sukshma')} for l,s,e,_ in sukshmas]),
                                            hide_index=True, use_container_width=True)
                                        if max_depth_sel >= 5:
                                            with st.expander("View Pranas", expanded=False):
                                                if sukshmas:
                                                    s_opt = [f"{p[0]} ({p[1].strftime('%Y-%m-%d %H:%M')} - {p[2].strftime('%Y-%m-%d %H:%M')})" for p in sukshmas]
                                                    sels = st.selectbox("Select Sukshma:", s_opt, key="sukshma_select")
                                                    pranas = sukshmas[s_opt.index(sels)][3]
                                                    st.dataframe(pd.DataFrame(
                                                        [{'Planet': l, 'Start': s.strftime('%Y-%m-%d %H:%M'),
                                                          'End': e.strftime('%Y-%m-%d %H:%M'),
                                                          'Duration': duration_str(e-s,'prana')} for l,s,e,_ in pranas]),
                                                        hide_index=True, use_container_width=True)
                                                    if max_depth_sel >= 6:
                                                        with st.expander("View Sub-Pranas", expanded=False):
                                                            if pranas:
                                                                p_opt = [f"{p[0]} ({p[1].strftime('%Y-%m-%d %H:%M')} - {p[2].strftime('%Y-%m-%d %H:%M')})" for p in pranas]
                                                                selp = st.selectbox("Select Prana:", p_opt, key="prana_select")
                                                                subp = pranas[p_opt.index(selp)][3]
                                                                st.dataframe(pd.DataFrame(
                                                                    [{'Planet': l, 'Start': s.strftime('%Y-%m-%d %H:%M'),
                                                                      'End': e.strftime('%Y-%m-%d %H:%M'),
                                                                      'Duration': duration_str(e-s,'sub_prana')} for l,s,e,_ in subp]),
                                                                    hide_index=True, use_container_width=True)
    st.info("Note: Periods are filtered from birth time; durations are approximate.")
    # =========================
    # SEPARATE INPUT: Current City → Live Micro-Periods
    # =========================
    st.subheader("Current City → Live Micro-Periods")
    cc1, cc2 = st.columns([2,1])
    with cc1:
        current_city_query = st.text_input(
            "Enter your CURRENT city (for local clock display)",
            placeholder="e.g., San Jose, CA or Chennai",
            key="current_city_input"
        )
    with cc2:
        depth_choice = st.selectbox(
            "Depth to inspect",
            ["Sukshma", "Prana", "Sub-Prana"],
            index=0
        )
    target_depth = _DEPTH_NAME_TO_INT[depth_choice]
    if st.button("Show current micro-periods", use_container_width=True):
        if not current_city_query or len(current_city_query) < 2:
            st.error("Please enter a valid current city.")
        else:
            try:
                cur_locs = geocode(current_city_query, exactly_one=False, limit=1)
                if not cur_locs:
                    st.error("Could not find that city. Try a more specific name.")
                else:
                    cur = cur_locs[0]
                    cur_lat, cur_lon = cur.latitude, cur.longitude
                    tz = tz_for_latlon(cur_lat, cur_lon)
                    # "Now" in that city's timezone, then convert to naive UTC for comparisons
                    now_local = datetime.now(tz)
                    now_utc_naive = now_local.astimezone(pytz.UTC).replace(tzinfo=None)
                    dp = cd['dasa_periods_filtered']
                    # Active stack Dasa→…→target_depth
                    active_path = find_active_path_to_depth(dp, now_utc_naive, target_depth)
                    flat_at_depth = collect_periods_at_depth(dp, target_depth)
                    # Locate current index at this depth
                    idx = None
                    for i, (_, s, e) in enumerate(flat_at_depth):
                        if s <= now_utc_naive < e:
                            idx = i
                            break
                    st.success(f"Detected time zone: {tz.zone} • Local now: {now_local.strftime('%Y-%m-%d %H:%M')}")
                    if active_path:
                        st.markdown("**Active stack (Dasa → … → " + depth_choice + ")**")
                        rows = []
                        labels = ["Dasa", "Bhukti", "Anthara", "Sukshma", "Prana", "Sub-Prana"]
                        level_key = depth_choice.lower().replace('-','_')
                        for j, (lord, s, e) in enumerate(active_path):
                            rows.append({
                                "Level": labels[j],
                                "Lord": lord,
                                "Start (local)": s.replace(tzinfo=pytz.UTC).astimezone(tz).strftime('%Y-%m-%d %H:%M'),
                                "End (local)": e.replace(tzinfo=pytz.UTC).astimezone(tz).strftime('%Y-%m-%d %H:%M'),
                                "Duration": duration_str(e - s, level=level_key if labels[j].lower()==depth_choice.lower() else 'dasa')
                            })
                        st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)
                    else:
                        st.info("Could not locate the active path at this depth. (Edge case near a boundary?)")
                    # Show current + next 5 at selected depth
                    if idx is not None:
                        nxt = flat_at_depth[idx: idx+6]
                        st.markdown(f"**Current and next {max(0, len(nxt)-1)} {depth_choice} periods (local time)**")
                        level_key = depth_choice.lower().replace('-','_')
                        tbl = []
                        for lord, s, e in nxt:
                            tbl.append({
                                "Lord": lord,
                                "Start (local)": s.replace(tzinfo=pytz.UTC).astimezone(tz).strftime('%Y-%m-%d %H:%M'),
                                "End (local)": e.replace(tzinfo=pytz.UTC).astimezone(tz).strftime('%Y-%m-%d %H:%M'),
                                "Duration": duration_str(e - s, level=level_key)
                            })
                        st.dataframe(pd.DataFrame(tbl), hide_index=True, use_container_width=True)
                    st.caption("Note: Micro-period boundaries are absolute (UTC). We display them in your city’s local clock.")
            except Exception as e:
                st.error(f"Could not compute local micro-periods: {e}")
else:
    st.info("Enter birth details above and click 'Generate Chart' to begin")
st.markdown("---")
st.caption("Sivapathy Astrology Data Generator")
