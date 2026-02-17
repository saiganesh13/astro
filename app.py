import streamlit as st
from datetime import datetime, timedelta
from math import sin, cos, tan, atan2, degrees, radians
from collections import defaultdict
import pandas as pd
import json
from geopy.geocoders import Nominatim
from geopy.extra.rate_limiter import RateLimiter
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from timezonefinder import TimezoneFinder
import pytz
import copy

# ---- Swiss Ephemeris / Astropy fallback ----
USE_SWISSEPH = False
try:
    import swisseph as swe
    import os as _os
    import urllib.request as _urlreq

    _EPHE_DIR = _os.path.join(_os.path.expanduser('~'), '.swisseph_ephe')
    _EPHE_BASE_URL = 'https://raw.githubusercontent.com/aloistr/swisseph/master/ephe/'
    _EPHE_FILES = ['sepl_18.se1', 'semo_18.se1', 'seas_18.se1']

    if not _os.path.isdir(_EPHE_DIR):
        _os.makedirs(_EPHE_DIR, exist_ok=True)

    for _ef in _EPHE_FILES:
        _ef_path = _os.path.join(_EPHE_DIR, _ef)
        if not _os.path.isfile(_ef_path):
            try:
                _urlreq.urlretrieve(_EPHE_BASE_URL + _ef, _ef_path)
            except Exception:
                pass  # will fallback to Moshier for missing files

    if any(_os.path.isfile(_os.path.join(_EPHE_DIR, f)) for f in _EPHE_FILES):
        swe.set_ephe_path(_EPHE_DIR)
    else:
        swe.set_ephe_path(None)  # fallback to built-in Moshier

    USE_SWISSEPH = True
except ImportError:
    from astropy.time import Time
    from astropy.coordinates import get_body, solar_system_ephemeris, GeocentricTrueEcliptic

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

# Sthana Bala Dict
sthana_bala_dict = {
    'Sun': [100,90,80,70,60,50,40,50,60,70,80,90],
    'Moon': [70,100,70,80,70,60,50,40,50,60,60,70], 
    'Jupiter': [60,60,70,100,90,60,75,60,80,40,50,80],
    'Venus': [60,70,60,50,40,35,80,50,60,80,70,100],
    'Mercury': [40,60,70,45,60,100,60,45,55,50,45,35],
    'Mars': [80,70,45,35,60,45,50,60,60,100,90,60],
    'Saturn': [35,50,60,70,80,60,100,90,50,60,80,50],
    'Rahu': [100]*12,
    'Ketu': [100]*12
}

# Status Mapping
status_data = {
    'Sun': {'Uchcham': 'Aries', 'Moolathirigonam': None, 'Aatchi': 'Leo', 'Neecham': 'Libra'},
    'Moon': {'Uchcham': 'Taurus', 'Moolathirigonam': None, 'Aatchi': 'Cancer', 'Neecham': 'Scorpio'},
    'Jupiter': {'Uchcham': 'Cancer', 'Moolathirigonam': 'Sagittarius', 'Aatchi': 'Pisces', 'Neecham': 'Capricorn'},
    'Venus': {'Uchcham': 'Pisces', 'Moolathirigonam': 'Libra', 'Aatchi': 'Taurus', 'Neecham': 'Virgo'},
    'Mercury': {'Uchcham': 'Virgo', 'Moolathirigonam': None, 'Aatchi': 'Gemini', 'Neecham': 'Pisces'},
    'Mars': {'Uchcham': 'Capricorn', 'Moolathirigonam': 'Aries', 'Aatchi': 'Scorpio', 'Neecham': 'Cancer'},
    'Saturn': {'Uchcham': 'Libra', 'Moolathirigonam': 'Aquarius', 'Aatchi': 'Capricorn', 'Neecham': 'Aries'}
}

capacity_dict = {
    'Saturn': 100, 'Mars': 100, 'Sun': 100, 'Jupiter': 100, 
    'Venus': 50, 'Mercury': 30, 'Moon': 100, 'Rahu': 100, 'Ketu': 50
}
good_capacity_dict = {
    'Saturn': 0, 'Mars': 25, 'Sun': 50, 'Jupiter': 100, 
    'Venus': 100, 'Mercury': 100, 'Rahu': 0, 'Ketu': 0
}
bad_capacity_dict = {
    'Saturn': 100, 'Mars': 75, 'Sun': 50, 'Jupiter': 0, 
    'Venus': 0, 'Mercury': 0, 'Rahu': 100, 'Ketu': 100
}

shukla_good = [100, 9, 16, 23, 30, 37, 44, 51, 58, 65, 72, 79, 86, 93, 100]
shukla_bad = [0] * 15
krishna_good = [93, 86, 79, 72, 65, 58, 51, 44, 37, 30, 23, 16, 9, 2, 0]
krishna_bad = [7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 84, 91, 98, 100]

shukla_tithi_names = ['Shukla Pratipada', 'Shukla Dwitiya', 'Shukla Tritiya', 'Shukla Chaturthi', 'Shukla Panchami', 'Shukla Shashti', 'Shukla Saptami', 'Shukla Ashtami', 'Shukla Navami', 'Shukla Dashami', 'Shukla Ekadashi', 'Shukla Dwadashi', 'Shukla Trayodashi', 'Shukla Chaturdashi', 'Purnima']
krishna_tithi_names = ['Krishna Pratipada', 'Krishna Dwitiya', 'Krishna Tritiya', 'Krishna Chaturthi', 'Krishna Panchami', 'Krishna Shashti', 'Krishna Saptami', 'Krishna Ashtami', 'Krishna Navami', 'Krishna Dashami', 'Krishna Ekadashi', 'Krishna Dwadashi', 'Krishna Trayodashi', 'Krishna Chaturdashi', 'Amavasya']

single_currency_planets = ['Venus', 'Jupiter', 'Mercury', 'Rahu', 'Ketu', 'Saturn']
bad_currency_planets = ['Saturn', 'Rahu', 'Ketu']
base_malefics = ['Saturn', 'Mars', 'Sun', 'Rahu']
malefic_planets = ['Saturn', 'Rahu', 'Ketu', 'Mars', 'Sun']

# MODIFICATION 3: Malefic Hierarchy for Navamsa Phase 1
navamsa_malefic_hierarchy = {'Rahu': 1, 'Saturn': 2, 'Sun': 3, 'Mars': 4, 'Ketu': 5}

mix_dict = {0:100,1:100,2:100,3:95,4:90,5:85,6:80,7:75,8:70,9:65,10:60,11:55,12:50,13:45,14:40,15:35,16:30,17:25,18:20,19:15,20:10,21:5,22:0}

# Maraivu Percentage by planet and house
maraivu_percentage = {
    'Sun':     {3: 25, 6: 75, 8: 100, 12: 50},
    'Moon':    {3: 25, 6: 75, 8: 100, 12: 50},
    'Venus':   {3: 100, 6: 25, 8: 100, 12: 0},
    'Mercury': {3: 25, 6: 50, 8: 100, 12: 50},
    'Jupiter': {3: 25, 6: 75, 8: 100, 12: 50},
    'Mars':    {3: 25, 6: 75, 8: 100, 12: 50},
    'Saturn':  {3: 25, 6: 75, 8: 100, 12: 50},
    'Rahu':    {3: 25, 6: 75, 8: 100, 12: 50},
    'Ketu':    {3: 25, 6: 75, 8: 100, 12: 50},
}

# MODIFICATION 1: Planet to ruled signs mapping
planet_ruled_signs = {
    'Sun': ['Leo'],
    'Moon': ['Cancer'],
    'Mars': ['Aries', 'Scorpio'],
    'Mercury': ['Gemini', 'Virgo'],
    'Jupiter': ['Sagittarius', 'Pisces'],
    'Venus': ['Taurus', 'Libra'],
    'Saturn': ['Capricorn', 'Aquarius']
}

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

def _datetime_to_jd(dt):
    """Convert a datetime to Julian Day, handling Julian/Gregorian calendar switch."""
    if USE_SWISSEPH:
        y, m, d = dt.year, dt.month, dt.day
        h = dt.hour + dt.minute / 60.0 + dt.second / 3600.0
        # Gregorian calendar from Oct 15 1582 onwards, Julian before
        if (y > 1582) or (y == 1582 and m > 10) or (y == 1582 and m == 10 and d >= 15):
            cal = swe.GREG_CAL
        else:
            cal = swe.JUL_CAL
        return swe.julday(y, m, d, h, cal)
    else:
        from astropy.time import Time as AstroTime
        return AstroTime(dt).jd

def compute_positions_swisseph(utc_dt, lat, lon):
    """Compute tropical planet longitudes + ascendant using Swiss Ephemeris."""
    # Set Lahiri ayanamsa mode BEFORE any ayanamsa query
    swe.set_sid_mode(swe.SIDM_LAHIRI)
    jd = _datetime_to_jd(utc_dt)

    # Planet IDs in swisseph
    planet_ids = {
        'sun': swe.SUN, 'moon': swe.MOON, 'mercury': swe.MERCURY,
        'venus': swe.VENUS, 'mars': swe.MARS, 'jupiter': swe.JUPITER,
        'saturn': swe.SATURN
    }

    lon_trop = {}
    for name, pid in planet_ids.items():
        result, _flag = swe.calc_ut(jd, pid)
        lon_trop[name] = result[0]  # tropical longitude

    # Rahu = True Node
    result, _flag = swe.calc_ut(jd, swe.TRUE_NODE)
    lon_trop['rahu'] = result[0]
    lon_trop['ketu'] = (result[0] + 180.0) % 360.0

    # Ascendant
    cusps, asmc = swe.houses(jd, lat, lon, b'P')  # Placidus
    asc_trop = asmc[0]

    # Use swisseph Lahiri ayanamsa for accuracy
    ayan_lahiri = swe.get_ayanamsa_ut(jd)

    return lon_trop, asc_trop, jd, ayan_lahiri

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

def get_navamsa_sign(lon):
    nav_lon = (lon * 9) % 360
    return get_sign(nav_lon)

def get_sign_lord(sign):
    sign_idx = sign_names.index(sign)
    return sign_lords[sign_idx]

def is_good_currency(c_key):
    if 'Bad' in c_key:
        return False
    if c_key in ['Jupiter', 'Venus', 'Mercury']:
        return True
    if c_key == 'Jupiter Poison':
        return True
    if 'Good' in c_key:
        return True
    return False

def is_sun_or_moon_currency(c_key):
    return 'Sun' in c_key or 'Moon' in c_key

def compute_chart(name, date_obj, time_str, lat, lon, tz_offset, max_depth):
    try:
        hour, minute = map(int, time_str.split(':'))
        if not (0<=hour<=23 and 0<=minute<=59): raise ValueError
    except:
        raise ValueError("Time must be in HH:MM format (24-hour)")
    
    local_dt = datetime.combine(date_obj, datetime.min.time().replace(hour=hour, minute=minute))
    utc_dt = local_dt - timedelta(hours=tz_offset)

    if USE_SWISSEPH:
        lon_trop, asc_trop, jd, ayan = compute_positions_swisseph(utc_dt, lat, lon)
        lon_sid = {p: get_sidereal_lon(v, ayan) for p, v in lon_trop.items()}
        lagna_sid = get_sidereal_lon(asc_trop, ayan)
    else:
        from astropy.time import Time
        from astropy.coordinates import get_body, solar_system_ephemeris, GeocentricTrueEcliptic
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
    
    if diff < 180:
        paksha = 'Shukla'
    else:
        paksha = 'Krishna'

    tithi_fraction = diff / 12
    tithi = int(tithi_fraction) + 1
    if tithi > 30: tithi = 30
    
    if paksha == 'Shukla':
        tithi_idx = tithi - 1
        if tithi_idx > 14: tithi_idx = 14
        moon_phase_name = shukla_tithi_names[tithi_idx]
    else:
        tithi_idx = tithi - 16
        if tithi_idx < 0: tithi_idx = 0
        if tithi_idx > 14: tithi_idx = 14
        moon_phase_name = krishna_tithi_names[tithi_idx]

    house_planets_rasi = defaultdict(list)
    positions = {**lon_sid, 'asc': lagna_sid}
    for p, L in positions.items():
        house_planets_rasi[get_house(L, lagna_sid)].append(p.capitalize() if p != 'asc' else 'Asc')

    planet_status_map = {}
    planet_sign_map = {}
    planet_house_map = {}
    for p in ['sun','moon','mars','mercury','jupiter','venus','saturn','rahu','ketu']:
        L = lon_sid[p]
        sign = get_sign(L)
        house = get_house(L, lagna_sid)
        planet_cap = p.capitalize()
        planet_sign_map[planet_cap] = sign
        planet_house_map[planet_cap] = house
        
        status = '-'
        if planet_cap in status_data:
            mapping = status_data[planet_cap]
            if sign == mapping['Uchcham']: status = 'Uchcham'
            elif sign == mapping['Neecham']: status = 'Neecham'
            elif sign == mapping['Moolathirigonam']: status = 'Moolathirigonam'
            elif sign == mapping['Aatchi']: status = 'Aatchi'
        planet_status_map[planet_cap] = status

    # MODIFICATION 1: Identify Uchcham planets and signs they rule
    uchcham_ruled_signs = set()
    for planet_cap, status in planet_status_map.items():
        if status == 'Uchcham' and planet_cap in planet_ruled_signs:
            for ruled_sign in planet_ruled_signs[planet_cap]:
                uchcham_ruled_signs.add(ruled_sign)
    
    parivardhana_map = {}
    planets_for_parivardhana = ['Sun', 'Moon', 'Mars', 'Mercury', 'Jupiter', 'Venus', 'Saturn']
    
    for planet_a in planets_for_parivardhana:
        sign_a = planet_sign_map[planet_a]
        lord_of_sign_a = get_sign_lord(sign_a)
        
        if lord_of_sign_a in planet_sign_map:
            sign_of_lord = planet_sign_map[lord_of_sign_a]
            lord_of_that_sign = get_sign_lord(sign_of_lord)
            
            if lord_of_that_sign == planet_a and planet_a != lord_of_sign_a:
                house_a = get_house(lon_sid[planet_a.lower()], lagna_sid)
                house_b = get_house(lon_sid[lord_of_sign_a.lower()], lagna_sid)
                parivardhana_map[planet_a] = f"{lord_of_sign_a} (H{house_a}-H{house_b})"
    
    rows = []
    planet_data = {}
    asc_deg = lagna_sid % 360; asc_sign = get_sign(asc_deg)
    a_nak, a_pada, a_ld, a_sl = get_nakshatra_details(asc_deg)
    dig_bala_asc = calculate_dig_bala('asc', asc_deg, lagna_sid)
    
    asc_nav_sign = get_navamsa_sign(asc_deg)
    asc_vargothuva = 'Yes' if asc_sign == asc_nav_sign else 'No'
    
    rows.append(['Asc', f"{asc_deg:.2f}", asc_sign, a_nak, a_pada, f"{a_ld}/{a_sl}", asc_vargothuva, '-',
                 f"{dig_bala_asc}%" if dig_bala_asc is not None else '', '', '', '', '', '', ''])
    
    moon_initial_good_val = 0.0
    
    for p in ['sun','moon','mars','mercury','jupiter','venus','saturn','rahu','ketu']:
        L = lon_sid[p]; sign = get_sign(L); nak, pada, ld, sl = get_nakshatra_details(L)
        dig_bala = calculate_dig_bala(p, L, lagna_sid)
        planet_cap = p.capitalize()
        sthana = sthana_bala_dict.get(planet_cap, [0]*12)[sign_names.index(sign)]
        
        nav_sign = get_navamsa_sign(L)
        vargothuva = 'Yes' if sign == nav_sign else 'No'
        parivardhana = parivardhana_map.get(planet_cap, '-')
        status = planet_status_map[planet_cap]
            
        capacity = capacity_dict.get(planet_cap, None)
        volume = (capacity * sthana / 100.0) if capacity is not None else 0.0
        
        # MODIFICATION 1: Apply +10% volume boost if planet is in a sign ruled by Uchcham planet
        if sign in uchcham_ruled_signs:
            volume = volume * 1.10
        
        moon_good_pct = 0
        moon_bad_pct = 0
        if planet_cap == 'Moon':
            if paksha == 'Shukla':
                good_pct = shukla_good[tithi_idx]
                bad_pct = shukla_bad[tithi_idx]
            else:
                good_pct = krishna_good[tithi_idx]
                bad_pct = krishna_bad[tithi_idx]
            moon_good_pct = good_pct
            moon_bad_pct = bad_pct
        else:
            good_pct = good_capacity_dict.get(planet_cap, 0)
            bad_pct = bad_capacity_dict.get(planet_cap, 0)

        good_val = volume * (good_pct / 100.0)
        bad_val = volume * (bad_pct / 100.0)
        
        if planet_cap == 'Moon':
            moon_initial_good_val = good_val
        
        total_debt = 0.0
        has_debt = False
        updated_status = '-'
        is_neechabhangam = False
        is_healthy_neecham_moon = False
        
        if status == 'Neecham':
            if planet_cap == 'Moon' and paksha == 'Shukla' and bad_val == 0:
                is_healthy_neecham_moon = True
            
            house_lord = get_sign_lord(sign)
            house_lord_status = planet_status_map.get(house_lord, '-')
            
            if house_lord_status in ['Uchcham', 'Moolathirigonam']:
                updated_status = 'Neechabhangam'
                is_neechabhangam = True
                
                nb_base_vol = capacity * 0.40
                neechabhangam_good_add = nb_base_vol * (good_pct / 100.0)
                neechabhangam_bad_add = nb_base_vol * (bad_pct / 100.0)
                
                good_val += neechabhangam_good_add
                bad_val += neechabhangam_bad_add
            else:
                current_house = planet_house_map[planet_cap]
                
                for other_planet in ['Sun', 'Moon', 'Mars', 'Mercury', 'Jupiter', 'Venus', 'Saturn', 'Rahu', 'Ketu']:
                    if other_planet != planet_cap:
                        other_house = planet_house_map.get(other_planet, -1)
                        if other_house == current_house:
                            other_status = planet_status_map.get(other_planet, '-')
                            if other_status in ['Uchcham', 'Moolathirigonam']:
                                updated_status = 'Neechabhangam'
                                is_neechabhangam = True
                                break
            
            if is_healthy_neecham_moon:
                if capacity is not None:
                    good_capacity = capacity * (good_pct / 100.0)
                    total_debt = -((1.2 * good_capacity) - good_val)
                    has_debt = True
            else:
                if capacity is not None:
                    total_debt = -((1.2 * capacity) - good_val)
                    has_debt = True
                    
        else:
            if bad_val > 0:
                total_debt = -bad_val
                has_debt = True

        if has_debt:
            debt_str = f"{total_debt:.2f}"
        else:
            debt_str = '-'
        
        currency_parts = []
        if planet_cap in single_currency_planets:
            total_val = good_val + bad_val
            if total_val > 0:
                if planet_cap in bad_currency_planets:
                    currency_parts.append(f"Bad {planet_cap}[{total_val:.2f}]")
                else:
                    currency_parts.append(f"{planet_cap}[{total_val:.2f}]")
        else:
            if good_val > 0:
                currency_parts.append(f"Good {planet_cap}[{good_val:.2f}]")
            if bad_val > 0:
                currency_parts.append(f"Bad {planet_cap}[{bad_val:.2f}]")

        default_currency_str = ", ".join(currency_parts)

        planet_data[planet_cap] = {
            'sthana': sthana, 'volume': volume, 'dig_bala': dig_bala, 'L': L, 
            'sign': sign, 'nak': nak, 'pada': pada, 'ld_sl': f"{ld}/{sl}", 
            'status': status, 'default_currency': default_currency_str,
            'debt': debt_str, 'vargothuva': vargothuva, 'updated_status': updated_status,
            'parivardhana': parivardhana,
            'good_inv': good_val, 'bad_inv': bad_val,
            'current_debt': total_debt,
            'final_inventory': defaultdict(float),
            'moon_phase': moon_phase_name if planet_cap == 'Moon' else None,
            'moon_bad_pct': moon_bad_pct if planet_cap == 'Moon' else 0,
            'moon_good_pct': moon_good_pct if planet_cap == 'Moon' else 0
        }
        
        if planet_cap in ['Jupiter', 'Venus', 'Mercury']:
            if good_val > 0: planet_data[planet_cap]['final_inventory'][planet_cap] = good_val
        elif planet_cap in ['Saturn', 'Rahu']:
            if bad_val > 0: planet_data[planet_cap]['final_inventory'][f"Bad {planet_cap}"] = bad_val
        elif planet_cap == 'Ketu':
            if good_val > 0: planet_data[planet_cap]['final_inventory'][f"Good {planet_cap}"] = good_val
            if bad_val > 0: planet_data[planet_cap]['final_inventory'][f"Bad {planet_cap}"] = bad_val
        else:
            if good_val > 0:
                key = f"Good {planet_cap}"
                if planet_cap == 'Moon': key = "Good Moon"
                planet_data[planet_cap]['final_inventory'][key] = good_val
            if bad_val > 0:
                key = f"Bad {planet_cap}"
                if planet_cap == 'Moon': key = "Bad Moon"
                planet_data[planet_cap]['final_inventory'][key] = bad_val

    # Mars in Leo: 100% Good Mars, no Bad Mars, no debt
    if planet_data['Mars']['sign'] == 'Leo':
        _mars_good = planet_data['Mars']['final_inventory'].get('Good Mars', 0.0)
        _mars_bad = planet_data['Mars']['final_inventory'].get('Bad Mars', 0.0)
        _mars_total = _mars_good + _mars_bad
        planet_data['Mars']['final_inventory']['Good Mars'] = _mars_total
        planet_data['Mars']['final_inventory'].pop('Bad Mars', None)
        planet_data['Mars']['current_debt'] = 0.0
        planet_data['Mars']['default_currency'] = f"Good Mars[{_mars_total:.2f}]"
        planet_data['Mars']['debt'] = "0.00"

    # Saturn in Taurus: switches from -100% malefic to -50 malefic and +50 benefic (Good Venus)
    if planet_data['Saturn']['sign'] == 'Taurus':
        _saturn_bad = planet_data['Saturn']['final_inventory'].get('Bad Saturn', 0.0)
        # Split: 50% stays as Bad Saturn, 50% becomes Venus (benefic)
        new_bad_saturn = _saturn_bad * 0.50
        new_good_venus = _saturn_bad * 0.50
        planet_data['Saturn']['final_inventory']['Bad Saturn'] = new_bad_saturn
        planet_data['Saturn']['final_inventory']['Venus'] = new_good_venus
        # Update debt to equal the new bad currency
        planet_data['Saturn']['current_debt'] = -new_bad_saturn if new_bad_saturn > 0 else 0.0
        # Update display strings
        saturn_parts = []
        if new_good_venus > 0: saturn_parts.append(f"Good Venus[{new_good_venus:.2f}]")
        if new_bad_saturn > 0: saturn_parts.append(f"Bad Saturn[{new_bad_saturn:.2f}]")
        planet_data['Saturn']['default_currency'] = ", ".join(saturn_parts)
        planet_data['Saturn']['debt'] = f"{planet_data['Saturn']['current_debt']:.2f}"

    # Rahu in Taurus: switches from -100% malefic to -50 malefic and +50 benefic (Venus currency)
    if planet_data['Rahu']['sign'] == 'Taurus':
        _rahu_bad = planet_data['Rahu']['final_inventory'].get('Bad Rahu', 0.0)
        # Split: 50% stays as Bad Rahu, 50% becomes Venus (benefic)
        new_bad_rahu = _rahu_bad * 0.50
        new_good_venus = _rahu_bad * 0.50
        planet_data['Rahu']['final_inventory']['Bad Rahu'] = new_bad_rahu
        planet_data['Rahu']['final_inventory']['Venus'] = new_good_venus
        # Update debt to equal the new bad currency
        planet_data['Rahu']['current_debt'] = -new_bad_rahu if new_bad_rahu > 0 else 0.0
        # Update display strings
        rahu_parts = []
        if new_good_venus > 0: rahu_parts.append(f"Good Venus[{new_good_venus:.2f}]")
        if new_bad_rahu > 0: rahu_parts.append(f"Bad Rahu[{new_bad_rahu:.2f}]")
        planet_data['Rahu']['default_currency'] = ", ".join(rahu_parts)
        planet_data['Rahu']['debt'] = f"{planet_data['Rahu']['current_debt']:.2f}"

    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        data = planet_data[p]
        rows.append([
            p, f"{data['L']:.2f}", data['sign'], data['nak'], data['pada'], data['ld_sl'], 
            data['vargothuva'],
            data['parivardhana'],
            f"{data['dig_bala']}%" if data['dig_bala'] is not None else '', 
            f"{data['sthana']}%", 
            data['status'],
            data['updated_status'],
            f"{data['volume']:.2f}", 
            data['default_currency'], 
            data['debt']
        ])
    
    # NAVAMSA EXCHANGE LOGIC (Phase 1) - MODIFICATION 3: Fixed Malefic Hierarchy
    nav_lagna = (lagna_sid * 9) % 360
    navamsa_data = {}
    navamsa_house_planets = defaultdict(list)
    nav_initial_default_debts = {}
    
    for p in ['sun','moon','mars','mercury','jupiter','venus','saturn','rahu','ketu']:
        L = lon_sid[p]
        nav_lon = (L * 9) % 360
        nav_house = (int(nav_lon/30) - int(nav_lagna/30)) % 12 + 1
        planet_cap = p.capitalize()
        navamsa_house_planets[nav_house].append(planet_cap)
        
        nav_volume = planet_data[planet_cap]['volume']
        
        if planet_cap == 'Moon':
            if paksha == 'Shukla':
                good_pct = shukla_good[tithi_idx]
                bad_pct = shukla_bad[tithi_idx]
            else:
                good_pct = krishna_good[tithi_idx]
                bad_pct = krishna_bad[tithi_idx]
            nav_moon_good_pct = good_pct
            nav_moon_bad_pct = bad_pct
        else:
            good_pct = good_capacity_dict.get(planet_cap, 0)
            bad_pct = bad_capacity_dict.get(planet_cap, 0)
            nav_moon_good_pct = 0
            nav_moon_bad_pct = 0
        
        nav_good_val = nav_volume * (good_pct / 100.0)
        nav_bad_val = nav_volume * (bad_pct / 100.0)
        
        nav_currency_parts = []
        if planet_cap in single_currency_planets:
            total_val = nav_good_val + nav_bad_val
            if total_val > 0:
                if planet_cap in bad_currency_planets:
                    nav_currency_parts.append(f"Bad {planet_cap}[{total_val:.2f}]")
                else:
                    nav_currency_parts.append(f"{planet_cap}[{total_val:.2f}]")
        else:
            if nav_good_val > 0:
                nav_currency_parts.append(f"Good {planet_cap}[{nav_good_val:.2f}]")
            if nav_bad_val > 0:
                nav_currency_parts.append(f"Bad {planet_cap}[{nav_bad_val:.2f}]")
        
        nav_default_currency_str = ", ".join(nav_currency_parts)
        
        nav_debt = 0.0
        if nav_bad_val > 0:
            nav_debt = -nav_bad_val
        
        nav_initial_default_debts[planet_cap] = nav_debt
        
        nav_inventory = defaultdict(float)
        if planet_cap in ['Jupiter', 'Venus', 'Mercury']:
            if nav_good_val > 0: nav_inventory[planet_cap] = nav_good_val
        elif planet_cap in ['Saturn', 'Rahu']:
            if nav_bad_val > 0: nav_inventory[f"Bad {planet_cap}"] = nav_bad_val
        elif planet_cap == 'Ketu':
            if nav_good_val > 0: nav_inventory[f"Good {planet_cap}"] = nav_good_val
            if nav_bad_val > 0: nav_inventory[f"Bad {planet_cap}"] = nav_bad_val
        else:
            if nav_good_val > 0:
                key = f"Good {planet_cap}"
                nav_inventory[key] = nav_good_val
            if nav_bad_val > 0:
                key = f"Bad {planet_cap}"
                nav_inventory[key] = nav_bad_val
        
        original_keys = set(nav_inventory.keys())
        
        navamsa_data[planet_cap] = {
            'nav_house': nav_house,
            'nav_volume': nav_volume,
            'nav_default_currency': nav_default_currency_str,
            'nav_good_val': nav_good_val,
            'nav_bad_val': nav_bad_val,
            'nav_current_debt': nav_debt,
            'nav_inventory': nav_inventory,
            'nav_moon_good_pct': nav_moon_good_pct,
            'nav_moon_bad_pct': nav_moon_bad_pct,
            'nav_original_keys': original_keys,
            'nav_gained_currencies': defaultdict(float),
            'nav_debt': 0.0
        }
    
    nav_debtor_rank = []
    nav_debtor_rank.append('Rahu')
    nav_debtor_rank.append('Sun')
    nav_debtor_rank.append('Saturn')
    
    nav_moon_p = navamsa_data['Moon']
    nav_is_waning = (paksha == 'Krishna')
    nav_bad_pct_moon = nav_moon_p['nav_moon_bad_pct']
    
    nav_moon_in_list = False
    
    if nav_is_waning and moon_phase_name == 'Amavasya':
        nav_debtor_rank.append('Moon')
        nav_moon_in_list = True
        
    nav_debtor_rank.append('Mars')
    
    if nav_is_waning and not nav_moon_in_list:
        if nav_bad_pct_moon > 25:
            nav_debtor_rank.append('Moon')
        else:
            nav_debtor_rank.append('Moon')
            
    nav_debtor_rank.append('Ketu')
    
    def get_nav_currency_rank_score(p_name, c_key):
        if p_name == 'Moon':
            phase = moon_phase_name
            is_shukla = (paksha == 'Shukla')
            idx = tithi_idx
            if phase == 'Purnima': return 1000
            if 'Good' in c_key or 'Moon' == c_key: 
                if is_shukla:
                    pct = shukla_good[idx]
                    return 500 + pct * 4 
                else:
                    pct = krishna_good[idx]
                    return 500 + pct * 4
        
        if c_key == 'Jupiter': return 990
        if c_key == 'Venus': return 980
        if c_key == 'Mercury': return 970
        if c_key == 'Good Mars': return 800 
        if c_key == 'Good Sun': return 700 
        if c_key == 'Good Ketu': return 700
        if p_name == 'Moon' and 'Bad' in c_key:
            pct = navamsa_data['Moon']['nav_moon_bad_pct']
            return 400 - (pct * 3)
        if c_key == 'Bad Mars': return 325 
        if c_key == 'Bad Sun': return 250 
        if c_key == 'Bad Saturn': return 100 
        if c_key == 'Bad Rahu': return 100
        if c_key == 'Bad Ketu': return 150
        return 0
    
    # MODIFICATION 3: Navamsa Exchange Cycle - Fixed Malefic Hierarchy
    nav_loop_active = True
    nav_cycle_limit = 200
    nav_cycles = 0
    
    while nav_loop_active and nav_cycles < nav_cycle_limit:
        nav_cycles += 1
        nav_something_happened = False
        
        for debtor in nav_debtor_rank:
            if navamsa_data[debtor]['nav_current_debt'] >= -0.001: continue
            
            debtor_house = navamsa_data[debtor]['nav_house']
            same_house_planets = navamsa_house_planets[debtor_house]
            
            debtor_is_malefic = debtor in malefic_planets
            debtor_malefic_rank = navamsa_malefic_hierarchy.get(debtor, 99)
            
            potential_targets = []
            for t_name in same_house_planets:
                if t_name == debtor: continue
                
                # MODIFICATION 3: Preserve Ketu logic - Ketu can only pull from Sun/Moon
                if debtor == 'Ketu' and t_name not in ['Sun', 'Moon']: continue
                
                # MODIFICATION 3: Allow malefic-on-malefic interactions based on hierarchy
                target_is_malefic = t_name in malefic_planets
                
                if debtor_is_malefic and target_is_malefic:
                    target_malefic_rank = navamsa_malefic_hierarchy.get(t_name, 99)
                    if debtor_malefic_rank >= target_malefic_rank:
                        continue
                else:
                    d_idx = nav_debtor_rank.index(debtor) if debtor in nav_debtor_rank else 99
                    t_idx = nav_debtor_rank.index(t_name) if t_name in nav_debtor_rank else 99
                    if d_idx > t_idx: continue 
                
                inv = navamsa_data[t_name]['nav_inventory']
                for key, val in inv.items():
                    if val > 0.001:
                        max_pull = navamsa_data[t_name]['nav_volume'] * 1.0
                        tracker_key = f"nav_pulled_from_{t_name}"
                        pulled = navamsa_data[debtor].get(tracker_key, 0.0)
                        
                        if pulled < max_pull:
                            score = get_nav_currency_rank_score(t_name, key)
                            is_good = is_good_currency(key)
                            potential_targets.append({
                                'planet': t_name, 'key': key, 'score': score, 'max_pull': max_pull,
                                'is_good': is_good
                            })

            if debtor_is_malefic:
                potential_targets.sort(key=lambda x: (-int(x['is_good']), -x['score']))
            else:
                potential_targets.sort(key=lambda x: -x['score'])
            
            good_available = any(t['is_good'] and navamsa_data[t['planet']]['nav_inventory'].get(t['key'], 0) > 0 for t in potential_targets)
            
            for tgt in potential_targets:
                if navamsa_data[debtor]['nav_current_debt'] >= -0.001: break
                
                if debtor_is_malefic and not tgt['is_good'] and good_available:
                    continue
                
                avail = navamsa_data[tgt['planet']]['nav_inventory'][tgt['key']]
                if avail <= 0: continue
                
                tracker_key = f"nav_pulled_from_{tgt['planet']}"
                pulled = navamsa_data[debtor].get(tracker_key, 0.0)
                cap_space = tgt['max_pull'] - pulled
                if cap_space <= 0: continue
                
                needed = abs(navamsa_data[debtor]['nav_current_debt'])
                take = min(1.0, avail, cap_space)
                
                if take > 0:
                    navamsa_data[tgt['planet']]['nav_inventory'][tgt['key']] -= take
                    navamsa_data[tgt['planet']]['nav_current_debt'] -= take
                    navamsa_data[tgt['planet']]['nav_debt'] -= take
                    
                    is_ketu_currency = (tgt['key'] == 'Bad Ketu' or tgt['key'] == 'Good Ketu')
                    is_sun_or_moon = (debtor in ['Sun', 'Moon'])
                    
                    if is_ketu_currency and not is_sun_or_moon:
                        navamsa_data[debtor]['nav_inventory']['Good Ketu'] += take
                        navamsa_data[debtor]['nav_gained_currencies']['Good Ketu'] += take
                        navamsa_data[debtor]['nav_current_debt'] += take
                        navamsa_data[debtor]['nav_debt'] += take
                    else:
                        navamsa_data[debtor]['nav_inventory'][tgt['key']] += take
                        navamsa_data[debtor]['nav_gained_currencies'][tgt['key']] += take
                        
                        is_bad_currency = 'Bad' in tgt['key'] or tgt['key'] in ['Amavasya', 'Bad Saturn', 'Bad Rahu']
                        if tgt['planet'] in ['Saturn', 'Rahu'] and 'Bad' in tgt['key']: is_bad_currency = True
                        if tgt['planet'] == 'Moon' and 'Bad' in tgt['key']: is_bad_currency = True
                        
                        if debtor == 'Ketu' and is_sun_or_moon_currency(tgt['key']):
                            navamsa_data[debtor]['nav_current_debt'] -= take
                            navamsa_data[debtor]['nav_debt'] -= take
                        elif is_bad_currency: 
                            navamsa_data[debtor]['nav_current_debt'] -= take
                            navamsa_data[debtor]['nav_debt'] -= take
                        else: 
                            navamsa_data[debtor]['nav_current_debt'] += take
                            navamsa_data[debtor]['nav_debt'] += take
                    
                    navamsa_data[debtor][tracker_key] = pulled + take
                    nav_something_happened = True
                    
                    good_available = any(t['is_good'] and navamsa_data[t['planet']]['nav_inventory'].get(t['key'], 0) > 0 for t in potential_targets)

        if not nav_something_happened: nav_loop_active = False
    
    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        gained = navamsa_data[p]['nav_gained_currencies']
        gained_parts = []
        for k, v in gained.items():
            if v > 0.001: gained_parts.append(f"{k}[{v:.2f}]")
        navamsa_data[p]['nav_gained_str'] = ", ".join(gained_parts) if gained_parts else "-"
        
        debt_val = navamsa_data[p]['nav_debt']
        if abs(debt_val) < 0.01:
            navamsa_data[p]['nav_debt_str'] = "0.00"
        else:
            navamsa_data[p]['nav_debt_str'] = f"{debt_val:.2f}"
    
    navamsa_rows = []
    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        nd = navamsa_data[p]
        navamsa_rows.append([
            p, 
            f"{nd['nav_volume']:.2f}", 
            nd['nav_default_currency'], 
            nd['nav_gained_str'],
            nd['nav_debt_str']
        ])
    
    df_navamsa_exchange = pd.DataFrame(navamsa_rows, columns=['Planet', 'Volume', 'Default Currency', 'Gained Currencies', 'Debt'])
    
    # NAVAMSA PHASE 2 LOGIC
    nav_moon_is_benefic_p2 = (paksha == 'Shukla') or (moon_phase_name == 'Purnima')
    nav_core_benefics_p2 = ['Jupiter', 'Venus', 'Mercury']
    if nav_moon_is_benefic_p2:
        nav_core_benefics_p2.append('Moon')
    
    navamsa_phase2_data = {}
    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        navamsa_phase2_data[p] = {
            'navp2_inventory': defaultdict(float),
            'navp2_current_debt': navamsa_data[p]['nav_current_debt'],
            'nav_volume': navamsa_data[p]['nav_volume'],
            'nav_house': navamsa_data[p]['nav_house'],
            'navp2_gained_currencies': defaultdict(float),
            'navp2_carried_over': defaultdict(float),
            'navp1_gained_currencies': defaultdict(float)
        }
        for k, v in navamsa_data[p]['nav_inventory'].items():
            navamsa_phase2_data[p]['navp2_inventory'][k] = v
            navamsa_phase2_data[p]['navp2_carried_over'][k] = v
        
        for k, v in navamsa_data[p]['nav_gained_currencies'].items():
            navamsa_phase2_data[p]['navp1_gained_currencies'][k] = v
    
    ketu_has_sun_moon = False
    for k in navamsa_phase2_data['Ketu']['navp2_inventory'].keys():
        if is_sun_or_moon_currency(k) and navamsa_phase2_data['Ketu']['navp2_inventory'][k] > 0.001:
            ketu_has_sun_moon = True
            break
    
    nav_bad_ketu_remaining = navamsa_phase2_data['Ketu']['navp2_inventory'].get('Bad Ketu', 0.0)
    if nav_bad_ketu_remaining > 0 and not ketu_has_sun_moon:
        navamsa_phase2_data['Ketu']['navp2_inventory']['Good Ketu'] = navamsa_phase2_data['Ketu']['navp2_inventory'].get('Good Ketu', 0.0) + nav_bad_ketu_remaining
        navamsa_phase2_data['Ketu']['navp2_carried_over']['Good Ketu'] = navamsa_phase2_data['Ketu']['navp2_carried_over'].get('Good Ketu', 0.0) + nav_bad_ketu_remaining
        navamsa_phase2_data['Ketu']['navp2_inventory']['Bad Ketu'] = 0.0
        navamsa_phase2_data['Ketu']['navp2_carried_over']['Bad Ketu'] = 0.0
    
    ketu_good_currency = navamsa_phase2_data['Ketu']['navp2_inventory'].get('Good Ketu', 0.0)
    if ketu_good_currency > 0 and not ketu_has_sun_moon:
        nav_core_benefics_p2.append('Ketu')
    
    nav_benefic_debt_pct = {}
    for p in nav_core_benefics_p2:
        volume = navamsa_phase2_data[p]['nav_volume']
        debt_p1 = abs(navamsa_phase2_data[p]['navp2_current_debt'])
        if volume > 0:
            debt_pct = (debt_p1 / volume) * 100
        else:
            debt_pct = 0.0
        nav_benefic_debt_pct[p] = debt_pct
    
    nav_sorted_benefics = sorted(nav_core_benefics_p2, key=lambda x: -nav_benefic_debt_pct[x])
    
    def get_navp2_currency_rank_score(p_name, c_key):
        if c_key == 'Jupiter': return 1000
        if c_key == 'Good Moon' or (p_name == 'Moon' and 'Good' in c_key): return 995
        if c_key == 'Venus': return 980
        if c_key == 'Mercury': return 970
        if c_key == 'Good Ketu': return 700
        return 0
    
    navp2_loop_active = True
    navp2_cycle_limit = 200
    navp2_cycles = 0
    
    while navp2_loop_active and navp2_cycles < navp2_cycle_limit:
        navp2_cycles += 1
        navp2_something_happened = False
        
        for puller in nav_sorted_benefics:
            if navamsa_phase2_data[puller]['navp2_current_debt'] >= -0.001: continue
            
            puller_debt_pct = nav_benefic_debt_pct[puller]
            puller_house = navamsa_phase2_data[puller]['nav_house']
            
            potential_targets = []
            for target in nav_core_benefics_p2:
                if target == puller: continue
                
                target_house = navamsa_phase2_data[target]['nav_house']
                if target_house != puller_house: continue
                
                if puller == 'Ketu' and target == 'Moon': continue
                if puller == 'Moon' and target == 'Ketu': continue
                
                target_debt_pct = nav_benefic_debt_pct[target]
                if target_debt_pct >= puller_debt_pct: continue
                
                inv = navamsa_phase2_data[target]['navp2_inventory']
                for key, val in inv.items():
                    if 'Bad' in key: continue
                    if val > 0.001:
                        max_pull = navamsa_phase2_data[target]['nav_volume'] * 1.0
                        tracker_key = f"navp2_pulled_from_{target}"
                        pulled = navamsa_phase2_data[puller].get(tracker_key, 0.0)
                        
                        if pulled < max_pull:
                            score = get_navp2_currency_rank_score(target, key)
                            potential_targets.append({
                                'planet': target, 'key': key, 'score': score, 'max_pull': max_pull
                            })
            
            potential_targets.sort(key=lambda x: -x['score'])
            
            for tgt in potential_targets:
                if navamsa_phase2_data[puller]['navp2_current_debt'] >= -0.001: break
                avail = navamsa_phase2_data[tgt['planet']]['navp2_inventory'][tgt['key']]
                if avail <= 0: continue
                
                tracker_key = f"navp2_pulled_from_{tgt['planet']}"
                pulled = navamsa_phase2_data[puller].get(tracker_key, 0.0)
                cap_space = tgt['max_pull'] - pulled
                if cap_space <= 0: continue
                
                take = min(1.0, avail, cap_space)
                
                if take > 0:
                    navamsa_phase2_data[tgt['planet']]['navp2_inventory'][tgt['key']] -= take
                    navamsa_phase2_data[tgt['planet']]['navp2_current_debt'] -= take
                    navamsa_phase2_data[puller]['navp2_inventory'][tgt['key']] += take
                    navamsa_phase2_data[puller]['navp2_gained_currencies'][tgt['key']] += take
                    navamsa_phase2_data[puller]['navp2_current_debt'] += take
                    navamsa_phase2_data[puller][tracker_key] = pulled + take
                    navp2_something_happened = True
                    
                    for ben in nav_core_benefics_p2:
                        vol = navamsa_phase2_data[ben]['nav_volume']
                        dbt = abs(navamsa_phase2_data[ben]['navp2_current_debt'])
                        if vol > 0:
                            nav_benefic_debt_pct[ben] = (dbt / vol) * 100
                        else:
                            nav_benefic_debt_pct[ben] = 0.0
        
        if not navp2_something_happened: navp2_loop_active = False
    
    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        inv = navamsa_phase2_data[p]['navp2_inventory']
        carried = navamsa_phase2_data[p]['navp2_carried_over']
        gained_p2 = navamsa_phase2_data[p]['navp2_gained_currencies']
        gained_p1 = navamsa_phase2_data[p]['navp1_gained_currencies']
        
        carried_parts = []
        own_keys = [p, f"Good {p}", f"Bad {p}"]
        if p == 'Moon': own_keys = ["Good Moon", "Bad Moon"]
        
        for k in carried.keys():
            current_val = inv.get(k, 0.0)
            gained_p2_val = gained_p2.get(k, 0.0)
            carried_val = current_val - gained_p2_val
            if carried_val > 0.001:
                carried_parts.append(f"{k}[{carried_val:.2f}]")
        
        navamsa_phase2_data[p]['currency_carried_str'] = ", ".join(carried_parts) if carried_parts else "-"
        
        combined_gained = defaultdict(float)
        for k, v in gained_p1.items():
            combined_gained[k] += v
        for k, v in gained_p2.items():
            combined_gained[k] += v
        
        navamsa_phase2_data[p]['combined_gained_p1_p2'] = combined_gained
        
        gained_parts = []
        for k, v in combined_gained.items():
            if v > 0.001: gained_parts.append(f"{k}[{v:.2f}]")
        navamsa_phase2_data[p]['currency_gained_str'] = ", ".join(gained_parts) if gained_parts else "-"
        
        initial_debt = nav_initial_default_debts.get(p, 0.0)
        corrected_debt = navamsa_phase2_data[p]['navp2_current_debt'] - initial_debt
        
        navamsa_phase2_data[p]['corrected_debt_p2'] = corrected_debt
        
        if abs(corrected_debt) < 0.01: 
            navamsa_phase2_data[p]['debt_navp2'] = "0.00"
        else: 
            navamsa_phase2_data[p]['debt_navp2'] = f"{corrected_debt:.2f}"
    
    navamsa_phase2_rows = []
    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        d_navp2 = navamsa_phase2_data[p]
        navamsa_phase2_rows.append([p, d_navp2['currency_carried_str'], d_navp2['currency_gained_str'], d_navp2['debt_navp2']])
    
    df_navamsa_phase2 = pd.DataFrame(navamsa_phase2_rows, columns=['Planet', 'Inventory Carried Over', 'Gained Currency', 'Debt [Nav Phase 2]'])
    
    # NAVAMSA PHASE 3 LOGIC - House Pot System
    navamsa_phase3_data = {}
    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        navamsa_phase3_data[p] = {
            'navp3_inventory': defaultdict(float),
            'navp3_current_debt': navamsa_phase2_data[p]['navp2_current_debt'],
            'nav_volume': navamsa_phase2_data[p]['nav_volume'],
            'nav_house': navamsa_phase2_data[p]['nav_house'],
            'navp3_debt': navamsa_phase2_data[p]['corrected_debt_p2'],
            'navp3_gained_currencies': defaultdict(float),
            'navp3_good_moon_gained': 0.0,
            'navp3_carried_over': defaultdict(float)
        }
        for k, v in navamsa_phase2_data[p]['navp2_inventory'].items():
            navamsa_phase3_data[p]['navp3_inventory'][k] = v
            navamsa_phase3_data[p]['navp3_carried_over'][k] = v
        
        for k, v in navamsa_phase2_data[p]['combined_gained_p1_p2'].items():
            navamsa_phase3_data[p]['navp3_gained_currencies'][k] = v
    
    house_pot = {}
    house_pot_allocations = {
        'Aries': 20, 'Taurus': 60, 'Gemini': 40, 'Cancer': moon_initial_good_val,
        'Leo': 30, 'Virgo': 60, 'Libra': 80, 'Scorpio': 20,
        'Sagittarius': 100, 'Capricorn': 10, 'Aquarius': 10, 'Pisces': 80
    }
    
    for h in range(1, 13):
        house_sign = get_sign((nav_lagna + (h - 1) * 30) % 360)
        house_pot[h] = house_pot_allocations.get(house_sign, 0)
    
    navp3_house_planets = defaultdict(list)
    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        h = navamsa_phase3_data[p]['nav_house']
        navp3_house_planets[h].append(p)
    
    navp3_cycle_limit = 200
    navp3_cycles = 0
    
    standard_malefics = ['Saturn', 'Rahu', 'Ketu', 'Mars', 'Sun']
    standard_benefics = ['Jupiter', 'Venus', 'Mercury']
    malefic_debtor_order = ['Rahu', 'Sun', 'Saturn', 'Mars', 'Ketu']
    
    while navp3_cycles < navp3_cycle_limit:
        navp3_cycles += 1
        navp3_something_happened = False
        
        for house_num in range(1, 13):
            if house_pot[house_num] <= 0.001:
                continue
            
            planets_in_house = navp3_house_planets[house_num]
            if not planets_in_house:
                continue
            
            moon_bad_currency = navamsa_phase3_data['Moon']['navp3_inventory'].get('Bad Moon', 0.0)
            moon_is_malefic = moon_bad_currency > 0.001
            
            house_malefics = []
            for p in planets_in_house:
                if p in standard_malefics:
                    house_malefics.append(p)
                elif p == 'Moon' and moon_is_malefic:
                    house_malefics.append(p)
            
            house_benefics = []
            for p in planets_in_house:
                if p in standard_benefics:
                    house_benefics.append(p)
                elif p == 'Moon' and not moon_is_malefic:
                    house_benefics.append(p)
            
            def get_malefic_rank(p):
                if p == 'Moon':
                    nav_vol = navamsa_phase3_data['Moon']['nav_volume']
                    if nav_vol > 0:
                        bad_pct = (moon_bad_currency / nav_vol) * 100
                    else:
                        bad_pct = 0
                    if bad_pct > 25:
                        return 2.5
                    else:
                        return 3.5
                else:
                    if p in malefic_debtor_order:
                        return malefic_debtor_order.index(p)
                    return 99
            
            house_malefics_sorted = sorted(house_malefics, key=get_malefic_rank)
            
            def get_benefic_debt_pct(p):
                vol = navamsa_phase3_data[p]['nav_volume']
                debt = abs(navamsa_phase3_data[p]['navp3_current_debt'])
                if vol > 0:
                    return (debt / vol) * 100
                return 0
            
            house_benefics_sorted = sorted(house_benefics, key=lambda p: -get_benefic_debt_pct(p))
            
            for malefic in house_malefics_sorted:
                if house_pot[house_num] <= 0.001:
                    break
                
                debt = navamsa_phase3_data[malefic]['navp3_current_debt']
                
                if debt < -0.001:
                    needed = abs(debt)
                    take = min(1.0, needed, house_pot[house_num])
                else:
                    take = min(1.0, house_pot[house_num])
                
                if take > 0.001:
                    house_pot[house_num] -= take
                    navamsa_phase3_data[malefic]['navp3_inventory']['Good Moon'] += take
                    navamsa_phase3_data[malefic]['navp3_current_debt'] += take
                    navamsa_phase3_data[malefic]['navp3_good_moon_gained'] += take
                    navamsa_phase3_data[malefic]['navp3_gained_currencies']['Good Moon'] += take
                    navamsa_phase3_data[malefic]['navp3_debt'] += take
                    navp3_something_happened = True
            
            all_malefics_cleared = True
            for malefic in house_malefics:
                if navamsa_phase3_data[malefic]['navp3_current_debt'] < -0.001:
                    all_malefics_cleared = False
                    break
            
            if all_malefics_cleared or len(house_malefics) == 0:
                for benefic in house_benefics_sorted:
                    if house_pot[house_num] <= 0.001:
                        break
                    
                    debt = navamsa_phase3_data[benefic]['navp3_current_debt']
                    
                    if debt < -0.001:
                        needed = abs(debt)
                        take = min(1.0, needed, house_pot[house_num])
                    else:
                        take = min(1.0, house_pot[house_num])
                    
                    if take > 0.001:
                        house_pot[house_num] -= take
                        navamsa_phase3_data[benefic]['navp3_inventory']['Good Moon'] += take
                        navamsa_phase3_data[benefic]['navp3_current_debt'] += take
                        navamsa_phase3_data[benefic]['navp3_good_moon_gained'] += take
                        navamsa_phase3_data[benefic]['navp3_gained_currencies']['Good Moon'] += take
                        navamsa_phase3_data[benefic]['navp3_debt'] += take
                        navp3_something_happened = True
        
        if not navp3_something_happened:
            break
    
    navamsa_phase3_rows = []
    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        inv = navamsa_phase3_data[p]['navp3_inventory']
        carried = navamsa_phase3_data[p]['navp3_carried_over']
        gained = navamsa_phase3_data[p]['navp3_gained_currencies']
        good_moon_p3 = navamsa_phase3_data[p]['navp3_good_moon_gained']
        
        carried_parts = []
        own_keys = [p, f"Good {p}", f"Bad {p}"]
        if p == 'Moon': own_keys = ["Good Moon", "Bad Moon"]
        
        for k in carried.keys():
            current_val = inv.get(k, 0.0)
            if k == 'Good Moon':
                carried_val = current_val - good_moon_p3
            else:
                carried_val = current_val
            if carried_val > 0.001:
                carried_parts.append(f"{k}[{carried_val:.2f}]")
        
        for k, v in inv.items():
            if k not in carried.keys() and k != 'Good Moon' and v > 0.001:
                carried_parts.append(f"{k}[{v:.2f}]")
        
        inventory_carried_str = ", ".join(carried_parts) if carried_parts else "-"
        
        gained_parts = []
        for k, v in gained.items():
            if v > 0.001: gained_parts.append(f"{k}[{v:.2f}]")
        gained_currencies_str = ", ".join(gained_parts) if gained_parts else "-"
        
        debt_val = navamsa_phase3_data[p]['navp3_debt']
        if abs(debt_val) < 0.01:
            debt_str = "0.00"
        else:
            debt_str = f"{debt_val:.2f}"
        
        navamsa_phase3_rows.append([p, inventory_carried_str, gained_currencies_str, debt_str])
    
    df_navamsa_phase3 = pd.DataFrame(navamsa_phase3_rows, columns=['Planet', 'Inventory Carried Over', 'Gained Currencies', 'Debt [Nav Phase 3]'])
    
    # PHASE 1 CURRENCY EXCHANGE LOGIC (Rasi Chart)
    debtor_rank = []
    debtor_rank.append('Rahu')
    debtor_rank.append('Sun')
    debtor_rank.append('Saturn')
    
    moon_p = planet_data['Moon']
    is_waning = (paksha == 'Krishna')
    bad_pct_moon = moon_p['moon_bad_pct']
    
    moon_in_list = False
    
    if is_waning and moon_p['moon_phase'] == 'Amavasya':
        debtor_rank.append('Moon')
        moon_in_list = True
        
    debtor_rank.append('Mars')
    
    if is_waning and not moon_in_list:
        if bad_pct_moon > 25:
            debtor_rank.append('Moon')
        else:
            debtor_rank.append('Moon')
            
    debtor_rank.append('Ketu')

    def get_currency_rank_score(p_name, c_key):
        if p_name == 'Moon':
            phase = planet_data['Moon']['moon_phase']
            is_shukla = (paksha == 'Shukla')
            idx = tithi_idx
            if phase == 'Purnima': return 1000
            if 'Good' in c_key or 'Moon' == c_key: 
                if is_shukla:
                    pct = shukla_good[idx]
                    return 500 + pct * 4 
                else:
                    pct = krishna_good[idx]
                    return 500 + pct * 4
        
        if c_key == 'Jupiter': return 990
        if c_key == 'Venus': return 980
        if c_key == 'Mercury': return 970
        if c_key == 'Good Mars': return 800 
        if c_key == 'Good Sun': return 700 
        if c_key == 'Good Ketu': return 700
        if p_name == 'Moon' and 'Bad' in c_key:
            pct = planet_data['Moon']['moon_bad_pct']
            return 400 - (pct * 3)
        if c_key == 'Bad Mars': return 325 
        if c_key == 'Bad Sun': return 250 
        if c_key == 'Bad Saturn': return 100 
        if c_key == 'Bad Rahu': return 100
        if c_key == 'Bad Ketu': return 150
        return 0

    loop_active = True
    cycle_limit = 200
    cycles = 0
    
    while loop_active and cycles < cycle_limit:
        cycles += 1
        something_happened = False
        
        for debtor in debtor_rank:
            if planet_data[debtor]['current_debt'] >= -0.001: continue
            
            debtor_is_malefic = debtor in malefic_planets
            
            potential_targets = []
            for t_name in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
                if t_name == debtor: continue
                d_idx = debtor_rank.index(debtor) if debtor in debtor_rank else 99
                t_idx = debtor_rank.index(t_name) if t_name in debtor_rank else 99
                if d_idx > t_idx: continue 
                if debtor == 'Ketu' and t_name not in ['Sun', 'Moon']: continue
                
                inv = planet_data[t_name]['final_inventory']
                for key, val in inv.items():
                    if val > 0.001:
                        L1 = planet_data[debtor]['L']
                        L2 = planet_data[t_name]['L']
                        diff = abs(L1 - L2)
                        if diff > 180: diff = 360 - diff
                        gap = int(diff)
                        if gap > 22: continue
                        
                        cap_pct = mix_dict.get(gap, 0)
                        max_pull = planet_data[t_name]['volume'] * (cap_pct / 100.0)
                        tracker_key = f"pulled_from_{t_name}"
                        pulled = planet_data[debtor].get(tracker_key, 0.0)
                        
                        if pulled < max_pull:
                            score = get_currency_rank_score(t_name, key)
                            is_good = is_good_currency(key)
                            potential_targets.append({
                                'planet': t_name, 'key': key, 'score': score, 'gap': gap, 'max_pull': max_pull,
                                'is_good': is_good
                            })

            if debtor_is_malefic:
                potential_targets.sort(key=lambda x: (-int(x['is_good']), -x['score'], x['gap']))
            else:
                potential_targets.sort(key=lambda x: (-x['score'], x['gap']))
            
            good_available = any(t['is_good'] and planet_data[t['planet']]['final_inventory'].get(t['key'], 0) > 0 for t in potential_targets)
            
            for tgt in potential_targets:
                if planet_data[debtor]['current_debt'] >= -0.001: break
                
                # Never pull your own bad currency from someone else
                _own_bad_key = f"Bad {debtor}" if debtor != 'Moon' else "Bad Moon"
                if tgt['key'] == _own_bad_key:
                    continue
                
                if debtor_is_malefic and not tgt['is_good'] and good_available:
                    continue
                
                avail = planet_data[tgt['planet']]['final_inventory'][tgt['key']]
                if avail <= 0: continue
                
                tracker_key = f"pulled_from_{tgt['planet']}"
                pulled = planet_data[debtor].get(tracker_key, 0.0)
                cap_space = tgt['max_pull'] - pulled
                if cap_space <= 0: continue
                
                needed = abs(planet_data[debtor]['current_debt'])
                take = min(1.0, avail, cap_space)
                
                if take > 0:
                    planet_data[tgt['planet']]['final_inventory'][tgt['key']] -= take
                    planet_data[tgt['planet']]['current_debt'] -= take
                    
                    is_ketu_currency = (tgt['key'] == 'Bad Ketu' or tgt['key'] == 'Good Ketu')
                    is_sun_or_moon = (debtor in ['Sun', 'Moon'])
                    
                    if is_ketu_currency and not is_sun_or_moon:
                        planet_data[debtor]['final_inventory']['Good Ketu'] += take
                        planet_data[debtor]['current_debt'] += take
                    else:
                        planet_data[debtor]['final_inventory'][tgt['key']] += take
                        
                        is_bad_currency_flag = 'Bad' in tgt['key'] or tgt['key'] in ['Amavasya', 'Bad Saturn', 'Bad Rahu']
                        if tgt['planet'] in ['Saturn', 'Rahu'] and 'Bad' in tgt['key']: is_bad_currency_flag = True
                        if tgt['planet'] == 'Moon' and 'Bad' in tgt['key']: is_bad_currency_flag = True
                        
                        if debtor == 'Ketu' and is_sun_or_moon_currency(tgt['key']):
                            planet_data[debtor]['current_debt'] -= take
                        elif is_bad_currency_flag: 
                            planet_data[debtor]['current_debt'] -= take 
                        else: 
                            planet_data[debtor]['current_debt'] += take
                    
                    planet_data[debtor][tracker_key] = pulled + take
                    something_happened = True
                    
                    # --- INFECTION PENALTY: Malefic-to-Malefic currency exchange ---
                    # If both Debtor (Receiver) and Giver (Target) are Malefic,
                    # inject the Debtor's Bad Currency into the Giver's inventory.
                    _tgt_name = tgt['planet']
                    _tgt_is_malefic = (_tgt_name in malefic_planets or
                                       (_tgt_name == 'Moon' and planet_data['Moon'].get('moon_bad_pct', 0) > 0))
                    if debtor_is_malefic and _tgt_is_malefic:
                        _infection_key = f"Bad {debtor}" if debtor != 'Moon' else "Bad Moon"
                        planet_data[_tgt_name]['final_inventory'][_infection_key] += take
                    
                    good_available = any(t['is_good'] and planet_data[t['planet']]['final_inventory'].get(t['key'], 0) > 0 for t in potential_targets)

        if not something_happened: loop_active = False

    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        inv = planet_data[p]['final_inventory']
        parts = []
        own_keys = [p, f"Good {p}", f"Bad {p}"]
        if p == 'Moon': own_keys = ["Good Moon", "Bad Moon"]
        for k in own_keys:
            if k in inv and inv[k] > 0.001: parts.append(f"{k}[{inv[k]:.2f}]")
        for k, v in inv.items():
            if k not in own_keys and v > 0.001: parts.append(f"{k}[{v:.2f}]")
        planet_data[p]['currency_p1'] = ", ".join(parts)
        d_val = planet_data[p]['current_debt']
        if abs(d_val) < 0.01: planet_data[p]['debt_p1'] = "0.00"
        else: planet_data[p]['debt_p1'] = f"{d_val:.2f}"

    phase1_rows = []
    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        d_p1 = planet_data[p]
        phase1_rows.append([p, d_p1['currency_p1'], d_p1['debt_p1']])
        
    df_phase1 = pd.DataFrame(phase1_rows, columns=['Planet', 'Currency [Phase 1]', 'Debt [Phase 1]'])
    
    # PHASE 2 CURRENCY EXCHANGE LOGIC (Rasi Chart)
    moon_is_benefic_p2 = (paksha == 'Shukla') or (moon_phase_name == 'Purnima')
    
    core_benefics_p2 = ['Jupiter', 'Venus', 'Mercury']
    if moon_is_benefic_p2:
        core_benefics_p2.append('Moon')
    
    malefics_p2 = ['Saturn', 'Rahu', 'Mars', 'Sun']
    
    phase2_data = {}
    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        phase2_data[p] = {
            'p2_inventory': defaultdict(float),
            'p2_current_debt': planet_data[p]['current_debt'],
            'volume': planet_data[p]['volume'],
            'L': planet_data[p]['L']
        }
        for k, v in planet_data[p]['final_inventory'].items():
            phase2_data[p]['p2_inventory'][k] = v
    
    rasi_ketu_has_sun_moon = False
    for k in phase2_data['Ketu']['p2_inventory'].keys():
        if is_sun_or_moon_currency(k) and phase2_data['Ketu']['p2_inventory'][k] > 0.001:
            rasi_ketu_has_sun_moon = True
            break
    
    bad_ketu_remaining = phase2_data['Ketu']['p2_inventory'].get('Bad Ketu', 0.0)
    if bad_ketu_remaining > 0 and not rasi_ketu_has_sun_moon:
        phase2_data['Ketu']['p2_inventory']['Good Ketu'] = phase2_data['Ketu']['p2_inventory'].get('Good Ketu', 0.0) + bad_ketu_remaining
        phase2_data['Ketu']['p2_inventory']['Bad Ketu'] = 0.0
        phase2_data['Ketu']['p2_current_debt'] += bad_ketu_remaining
    
    if not rasi_ketu_has_sun_moon:
        core_benefics_p2.append('Ketu')
    
    benefic_debt_pct = {}
    for p in core_benefics_p2:
        volume = phase2_data[p]['volume']
        debt_p1 = abs(phase2_data[p]['p2_current_debt'])
        if volume > 0:
            debt_pct = (debt_p1 / volume) * 100
        else:
            debt_pct = 0.0
        benefic_debt_pct[p] = debt_pct
    
    sorted_benefics = sorted(core_benefics_p2, key=lambda x: -benefic_debt_pct[x])
    
    def get_p2_currency_rank_score(p_name, c_key):
        if c_key == 'Jupiter': return 1000
        if c_key == 'Good Moon' or (p_name == 'Moon' and 'Good' in c_key): return 995
        if c_key == 'Venus': return 980
        if c_key == 'Mercury': return 970
        if c_key == 'Good Ketu': return 700
        return 0
    
    p2_loop_active = True
    p2_cycle_limit = 200
    p2_cycles = 0
    
    while p2_loop_active and p2_cycles < p2_cycle_limit:
        p2_cycles += 1
        p2_something_happened = False
        
        for puller in sorted_benefics:
            if phase2_data[puller]['p2_current_debt'] >= -0.001: continue
            
            puller_debt_pct = benefic_debt_pct[puller]
            
            potential_targets = []
            for target in core_benefics_p2:
                if target == puller: continue
                
                if puller == 'Ketu' and target == 'Moon': continue
                if puller == 'Moon' and target == 'Ketu': continue
                
                target_debt_pct = benefic_debt_pct[target]
                if target_debt_pct >= puller_debt_pct: continue
                
                L1 = phase2_data[puller]['L']
                L2 = phase2_data[target]['L']
                diff = abs(L1 - L2)
                if diff > 180: diff = 360 - diff
                gap = int(diff)
                if gap > 22: continue
                
                inv = phase2_data[target]['p2_inventory']
                for key, val in inv.items():
                    if 'Bad' in key: continue
                    if val > 0.001:
                        cap_pct = mix_dict.get(gap, 0)
                        max_pull = phase2_data[target]['volume'] * (cap_pct / 100.0)
                        tracker_key = f"p2_pulled_from_{target}"
                        pulled = phase2_data[puller].get(tracker_key, 0.0)
                        
                        if pulled < max_pull:
                            score = get_p2_currency_rank_score(target, key)
                            potential_targets.append({
                                'planet': target, 'key': key, 'score': score, 'gap': gap, 'max_pull': max_pull
                            })
            
            potential_targets.sort(key=lambda x: (-x['score'], x['gap']))
            
            for tgt in potential_targets:
                if phase2_data[puller]['p2_current_debt'] >= -0.001: break
                avail = phase2_data[tgt['planet']]['p2_inventory'][tgt['key']]
                if avail <= 0: continue
                
                tracker_key = f"p2_pulled_from_{tgt['planet']}"
                pulled = phase2_data[puller].get(tracker_key, 0.0)
                cap_space = tgt['max_pull'] - pulled
                if cap_space <= 0: continue
                
                take = min(1.0, avail, cap_space)
                
                if take > 0:
                    phase2_data[tgt['planet']]['p2_inventory'][tgt['key']] -= take
                    phase2_data[tgt['planet']]['p2_current_debt'] -= take
                    phase2_data[puller]['p2_inventory'][tgt['key']] += take
                    phase2_data[puller]['p2_current_debt'] += take
                    phase2_data[puller][tracker_key] = pulled + take
                    p2_something_happened = True
                    
                    for ben in core_benefics_p2:
                        vol = phase2_data[ben]['volume']
                        dbt = abs(phase2_data[ben]['p2_current_debt'])
                        if vol > 0:
                            benefic_debt_pct[ben] = (dbt / vol) * 100
                        else:
                            benefic_debt_pct[ben] = 0.0
        
        if not p2_something_happened: p2_loop_active = False
    
    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        inv = phase2_data[p]['p2_inventory']
        parts = []
        own_keys = [p, f"Good {p}", f"Bad {p}"]
        if p == 'Moon': own_keys = ["Good Moon", "Bad Moon"]
        for k in own_keys:
            if k in inv and inv[k] > 0.001: parts.append(f"{k}[{inv[k]:.2f}]")
        for k, v in inv.items():
            if k not in own_keys and v > 0.001: parts.append(f"{k}[{v:.2f}]")
        phase2_data[p]['currency_p2'] = ", ".join(parts) if parts else "-"
        d_val = phase2_data[p]['p2_current_debt']
        if abs(d_val) < 0.01: phase2_data[p]['debt_p2'] = "0.00"
        else: phase2_data[p]['debt_p2'] = f"{d_val:.2f}"
    
    phase2_rows = []
    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        d_p2 = phase2_data[p]
        phase2_rows.append([p, d_p2['currency_p2'], d_p2['debt_p2']])
    
    df_phase2 = pd.DataFrame(phase2_rows, columns=['Planet', 'Currency [Phase 2]', 'Debt [Phase 2]'])
    
    # GLOBAL RESERVE INITIALIZATION
    house_reserves = defaultdict(lambda: defaultdict(float))
    
    # PHASE 3 CURRENCY EXCHANGE LOGIC (Rasi Chart) - 11th House Pot System
    phase3_data = {}
    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        phase3_data[p] = {
            'p3_inventory': defaultdict(float),
            'p3_current_debt': phase2_data[p]['p2_current_debt'],
            'volume': phase2_data[p]['volume'],
            'L': phase2_data[p]['L'],
            'rasi_house': planet_house_map[p]
        }
        for k, v in phase2_data[p]['p2_inventory'].items():
            phase3_data[p]['p3_inventory'][k] = v
    
    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        navp3_gained = navamsa_phase3_data[p]['navp3_gained_currencies']
        
        for k, v in navp3_gained.items():
            if v > 0.001:
                add_amount = v * 0.10
                phase3_data[p]['p3_inventory'][k] += add_amount
        
        navp3_debt = navamsa_phase3_data[p]['navp3_debt']
        phase3_data[p]['p3_current_debt'] += navp3_debt * 0.10
    
    planets_in_house_11 = []
    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        if phase3_data[p]['rasi_house'] == 11:
            planets_in_house_11.append(p)
    
    house_11_pot = 50.0
    house_11_sign = get_sign((lagna_sid + (11 - 1) * 30) % 360)
    
    if planets_in_house_11:
        # Collect planets that have debt (negative p3_current_debt)
        debtor_planets = []
        for p in planets_in_house_11:
            if phase3_data[p]['p3_current_debt'] < -0.001:
                debtor_planets.append(p)
        
        if debtor_planets:
            # Each planet's absolute debt and gift cap (50% of its debt)
            planet_debts = {}
            planet_caps = {}
            for p in debtor_planets:
                abs_debt = abs(phase3_data[p]['p3_current_debt'])
                planet_debts[p] = abs_debt
                planet_caps[p] = abs_debt * 0.50  # 50% of debt is the max gift cap
            
            # Distribute the pot proportionally based on debt, respecting 50% cap
            remaining_pot = house_11_pot
            remaining_planets = list(debtor_planets)
            gifts = {p: 0.0 for p in debtor_planets}
            
            while remaining_pot > 0.001 and remaining_planets:
                remaining_total_debt = sum(planet_debts[p] for p in remaining_planets)
                if remaining_total_debt < 0.001:
                    break
                
                # Calculate proportional shares for this round
                shares = {}
                for p in remaining_planets:
                    shares[p] = (planet_debts[p] / remaining_total_debt) * remaining_pot
                
                # Find planets whose proportional share would exceed their 50% cap
                capped = []
                for p in remaining_planets:
                    remaining_cap = planet_caps[p] - gifts[p]
                    if shares[p] >= remaining_cap - 0.001:
                        capped.append(p)
                
                if not capped:
                    # No planet hit its cap - distribute all proportional shares
                    for p in remaining_planets:
                        gifts[p] += shares[p]
                    remaining_pot = 0.0
                    break
                else:
                    # Give capped planets their cap amount, then redistribute the rest
                    for p in capped:
                        actual_gift = max(0, planet_caps[p] - gifts[p])
                        gifts[p] += actual_gift
                        remaining_pot -= actual_gift
                    remaining_planets = [p for p in remaining_planets if p not in capped]
            
            # Apply the computed gifts to each planet
            for p in debtor_planets:
                if gifts[p] > 0.001:
                    phase3_data[p]['p3_inventory']['Good Moon'] += gifts[p]
                    phase3_data[p]['p3_current_debt'] += gifts[p]
            
            house_11_pot = max(0, remaining_pot)
    
    if house_11_pot > 0.001:
        house_reserves[house_11_sign]['Good Moon'] += house_11_pot
    
    phase3_rows = []
    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        inv = phase3_data[p]['p3_inventory']
        parts = []
        own_keys = [p, f"Good {p}", f"Bad {p}"]
        if p == 'Moon': own_keys = ["Good Moon", "Bad Moon"]
        for k in own_keys:
            if k in inv and inv[k] > 0.001: parts.append(f"{k}[{inv[k]:.2f}]")
        for k, v in inv.items():
            if k not in own_keys and v > 0.001: parts.append(f"{k}[{v:.2f}]")
        phase3_data[p]['currency_p3'] = ", ".join(parts) if parts else "-"
        d_val = phase3_data[p]['p3_current_debt']
        if abs(d_val) < 0.01: phase3_data[p]['debt_p3'] = "0.00"
        else: phase3_data[p]['debt_p3'] = f"{d_val:.2f}"
    
    phase3_rows = []
    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        d_p3 = phase3_data[p]
        phase3_rows.append([p, d_p3['currency_p3'], d_p3['debt_p3']])
    
    df_phase3 = pd.DataFrame(phase3_rows, columns=['Planet', 'Currency [Phase 3]', 'Debt [Phase 3]'])
    
    # PHASE 4 CURRENCY EXCHANGE LOGIC (Rasi Chart) - Gift Pots System
    phase4_data = {}
    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        phase4_data[p] = {
            'p4_inventory': defaultdict(float),
            'p4_current_debt': phase3_data[p]['p3_current_debt'],
            'volume': phase3_data[p]['volume'],
            'L': phase3_data[p]['L'],
            'rasi_house': phase3_data[p]['rasi_house'],
            'sign': planet_sign_map[p]
        }
        for k, v in phase3_data[p]['p3_inventory'].items():
            phase4_data[p]['p4_inventory'][k] = v
    
    gift_pot_config = {
        'Sagittarius': ('Jupiter', 100),
        'Pisces': ('Jupiter', 80),
        'Libra': ('Venus', 80),
        'Taurus': ('Venus', 60)
    }
    
    pot_inventory = {}
    pot_currency_type = {}
    
    for sign_name, (gifter, multiplier) in gift_pot_config.items():
        gifter_sthana = planet_data[gifter]['sthana']
        pot_value = multiplier * (gifter_sthana / 100.0)
        pot_inventory[sign_name] = pot_value
        pot_currency_type[sign_name] = gifter
    
    p4_standard_malefics = ['Saturn', 'Rahu', 'Ketu', 'Mars', 'Sun']
    p4_standard_benefics = ['Jupiter', 'Venus', 'Mercury']
    p4_malefic_debtor_order = ['Rahu', 'Sun', 'Saturn', 'Mars', 'Ketu']
    
    for target_sign in ['Sagittarius', 'Pisces', 'Libra', 'Taurus']:
        sign_pot = pot_inventory[target_sign]
        currency_type = pot_currency_type[target_sign]
        
        if sign_pot <= 0.001:
            continue
        
        planets_in_sign = []
        for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
            if phase4_data[p]['sign'] == target_sign:
                planets_in_sign.append(p)
        
        if not planets_in_sign:
            house_reserves[target_sign][currency_type] += sign_pot
            pot_inventory[target_sign] = 0.0
            continue
        
        p4_cycle_limit = 200
        p4_cycles = 0
        
        while p4_cycles < p4_cycle_limit:
            p4_cycles += 1
            p4_something_happened = False
            
            if sign_pot <= 0.001:
                break
            
            moon_bad_currency_p4 = phase4_data['Moon']['p4_inventory'].get('Bad Moon', 0.0)
            moon_is_malefic_p4 = moon_bad_currency_p4 > 0.001
            
            sign_malefics = []
            for p in planets_in_sign:
                if p in p4_standard_malefics:
                    sign_malefics.append(p)
                elif p == 'Moon' and moon_is_malefic_p4:
                    sign_malefics.append(p)
            
            sign_benefics = []
            for p in planets_in_sign:
                if p in p4_standard_benefics:
                    sign_benefics.append(p)
                elif p == 'Moon' and not moon_is_malefic_p4:
                    sign_benefics.append(p)
            
            def get_p4_malefic_rank(p):
                if p == 'Moon':
                    vol = phase4_data['Moon']['volume']
                    if vol > 0:
                        bad_pct = (moon_bad_currency_p4 / vol) * 100
                    else:
                        bad_pct = 0
                    if bad_pct > 25:
                        return 2.5
                    else:
                        return 3.5
                else:
                    if p in p4_malefic_debtor_order:
                        return p4_malefic_debtor_order.index(p)
                    return 99
            
            sign_malefics_sorted = sorted(sign_malefics, key=get_p4_malefic_rank)
            
            def get_p4_benefic_debt_pct(p):
                vol = phase4_data[p]['volume']
                debt = abs(phase4_data[p]['p4_current_debt'])
                if vol > 0:
                    return (debt / vol) * 100
                return 0
            
            sign_benefics_sorted = sorted(sign_benefics, key=lambda p: -get_p4_benefic_debt_pct(p))
            
            for malefic in sign_malefics_sorted:
                if sign_pot <= 0.001:
                    break
                
                debt = phase4_data[malefic]['p4_current_debt']
                
                if debt < -0.001:
                    needed = abs(debt)
                    take = min(1.0, needed, sign_pot)
                    
                    if take > 0.001:
                        sign_pot -= take
                        phase4_data[malefic]['p4_inventory'][currency_type] += take
                        phase4_data[malefic]['p4_current_debt'] += take
                        p4_something_happened = True
            
            all_malefics_cleared = True
            for malefic in sign_malefics:
                if phase4_data[malefic]['p4_current_debt'] < -0.001:
                    all_malefics_cleared = False
                    break
            
            if all_malefics_cleared or len(sign_malefics) == 0:
                for benefic in sign_benefics_sorted:
                    if sign_pot <= 0.001:
                        break
                    
                    debt = phase4_data[benefic]['p4_current_debt']
                    
                    if debt < -0.001:
                        needed = abs(debt)
                        take = min(1.0, needed, sign_pot)
                        
                        if take > 0.001:
                            sign_pot -= take
                            phase4_data[benefic]['p4_inventory'][currency_type] += take
                            phase4_data[benefic]['p4_current_debt'] += take
                            p4_something_happened = True
            
            if not p4_something_happened:
                break
        
        pot_inventory[target_sign] = sign_pot
        
        if sign_pot > 0.001:
            house_reserves[target_sign][currency_type] += sign_pot
    
    phase4_rows = []
    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        inv = phase4_data[p]['p4_inventory']
        parts = []
        own_keys = [p, f"Good {p}", f"Bad {p}"]
        if p == 'Moon': own_keys = ["Good Moon", "Bad Moon"]
        for k in own_keys:
            if k in inv and inv[k] > 0.001: parts.append(f"{k}[{inv[k]:.2f}]")
        for k, v in inv.items():
            if k not in own_keys and v > 0.001: parts.append(f"{k}[{v:.2f}]")
        phase4_data[p]['currency_p4'] = ", ".join(parts) if parts else "-"
        d_val = phase4_data[p]['p4_current_debt']
        if abs(d_val) < 0.01: phase4_data[p]['debt_p4'] = "0.00"
        else: phase4_data[p]['debt_p4'] = f"{d_val:.2f}"
    
    phase4_rows = []
    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        d_p4 = phase4_data[p]
        phase4_rows.append([p, d_p4['currency_p4'], d_p4['debt_p4']])
    
    df_phase4 = pd.DataFrame(phase4_rows, columns=['Planet', 'Currency [Phase 4]', 'Debt [Phase 4]'])
    
    # PHASE 5 CURRENCY EXCHANGE LOGIC - Virtual Aspect Clones
    # MODIFICATION 2: Reordered Steps
    # Step 1: Virtual Malefic Clones Pull (Active Pulling) - currency marked as "Wasted"
    # Step 2: Real Malefics Pull (from Clone's original inventory only)
    # Step 3: Real Benefics Pull (from Clone's original inventory only)
    
    phase5_data = {}
    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        phase5_data[p] = {
            'p5_inventory': defaultdict(float),
            'p5_current_debt': phase4_data[p]['p4_current_debt'],
            'volume': phase4_data[p]['volume'],
            'L': phase4_data[p]['L'],
            'rasi_house': phase4_data[p]['rasi_house'],
            'sign': phase4_data[p]['sign'],
            'bad_inv': 0.0
        }
        for k, v in phase4_data[p]['p4_inventory'].items():
            phase5_data[p]['p5_inventory'][k] = v
            if 'Bad' in k:
                phase5_data[p]['bad_inv'] += v
    
    P5_STANDARD_MALEFICS = ['Saturn', 'Mars', 'Sun', 'Rahu', 'Ketu']
    P5_STANDARD_BENEFICS = ['Jupiter', 'Venus', 'Mercury']
    
    ASPECT_RULES = {
        'Saturn': {3: 0.25, 7: 1.0, 10: 0.75},
        'Mars': {4: 0.40, 7: 1.0, 8: 0.25},
        'Jupiter': {5: 1.0, 7: 1.0, 9: 1.0},
        'Venus': {7: 1.0},
        'Mercury': {7: 1.0},
        'Moon': {4: 0.25, 6: 0.50, 7: 1.0, 8: 0.50, 10: 0.25}
    }
    
    P5_MALEFIC_DEBTOR_RANK = ['Rahu', 'Sun', 'Saturn', 'Mars', 'Ketu']
    
    leftover_aspects = []
    all_leftover_clones = []
    all_initial_clones = []
    
    def get_p5_currency_rank_score(c_key):
        if c_key == 'Jupiter': return 990
        if c_key == 'Jupiter Poison': return 990
        if c_key == 'Venus': return 980
        if c_key == 'Mercury': return 970
        if c_key == 'Good Moon': return 950
        if c_key == 'Good Saturn': return 780
        if c_key == 'Good Mars': return 770
        if c_key == 'Good Sun': return 760
        if c_key == 'Good Ketu': return 700
        if c_key == 'Bad Moon': return 300
        if c_key == 'Bad Mars': return 250
        if c_key == 'Bad Sun': return 200
        if c_key == 'Bad Saturn': return 100
        if c_key == 'Bad Rahu': return 100
        if c_key == 'Bad Ketu': return 150
        return 0
    
    def is_moon_malefic_p5():
        return phase5_data['Moon']['bad_inv'] > 0.001
    
    PLANET_SEQUENCE = ['Saturn', 'Mars', 'Jupiter', 'Venus', 'Mercury', 'Moon']
    _jp_poison_notes = []  # Jupiter Poison diagnostic notes (initialized before loop)
    
    for current_planet in PLANET_SEQUENCE:
        if current_planet not in ASPECT_RULES:
            continue
        
        aspect_offsets = ASPECT_RULES[current_planet]
        parent_data = phase5_data[current_planet]
        parent_L = parent_data['L']
        parent_debt = parent_data['p5_current_debt']

        # Local scaling for Malefic Neecham/Neechabhangam planets (applied to clones only)
        _cp_is_malefic = (
            current_planet in ('Saturn', 'Mars', 'Sun', 'Rahu', 'Ketu')
            or (current_planet == 'Moon' and phase5_data['Moon']['bad_inv'] > 0.001)
        )
        _cp_status = planet_data[current_planet].get('updated_status') or planet_data[current_planet].get('status', '')
        if _cp_is_malefic and _cp_status in ('Neecham', 'Neechabhangam', 'Neechabhanga Raja Yoga'):
            scaling_factor = planet_data[current_planet]['sthana'] / 120.0
        else:
            scaling_factor = 1.0
        
        # Part A: Clone Creation
        clones = []
        
        for offset, aspect_pct in aspect_offsets.items():
            clone_L = (parent_L + (offset - 1) * 30) % 360
            parent_inv = parent_data['p5_inventory']
            
            clone_inventory = defaultdict(float)
            clone_debt = 0.0
            clone_type = 'Passive'
            
            if current_planet == 'Saturn':
                good_sum = 0.0
                for k, v in parent_inv.items():
                    if is_good_currency(k) and v > 0.001:
                        good_sum += v / 2.0
                
                clone_value = aspect_pct * good_sum * scaling_factor
                if clone_value > 0.001:
                    clone_inventory['Good Saturn'] = clone_value
                
                clone_debt = parent_debt * aspect_pct * scaling_factor
                clone_type = 'Active'
                
            elif current_planet == 'Mars':
                good_sum = 0.0
                for k, v in parent_inv.items():
                    if is_good_currency(k) and v > 0.001:
                        if k == 'Good Mars':
                            good_sum += v
                        else:
                            good_sum += v / 2.0
                
                clone_value = aspect_pct * good_sum * scaling_factor
                if clone_value > 0.001:
                    clone_inventory['Good Mars'] = clone_value
                
                clone_debt = parent_debt * aspect_pct * scaling_factor
                clone_type = 'Active'
                
            elif current_planet in ['Jupiter', 'Venus', 'Mercury']:
                good_sum = 0.0
                for k, v in parent_inv.items():
                    if is_good_currency(k) and v > 0.001:
                        good_sum += v
                
                clone_value = aspect_pct * good_sum * scaling_factor
                if clone_value > 0.001:
                    clone_inventory[current_planet] = clone_value
                
                clone_debt = 0.0
                clone_type = 'Passive'
                
            elif current_planet == 'Moon':
                good_moon_val = parent_inv.get('Good Moon', 0.0)
                other_good_sum = 0.0
                for k, v in parent_inv.items():
                    if is_good_currency(k) and k != 'Good Moon' and v > 0.001:
                        other_good_sum += v
                
                total_value = good_moon_val + (other_good_sum / 2.0)
                clone_value = aspect_pct * total_value * scaling_factor
                if clone_value > 0.001:
                    clone_inventory['Good Moon'] = clone_value
                
                clone_debt = 0.0
                clone_type = 'Passive'
            
            # MODIFICATION 2: Store original inventory and wasted inventory separately
            original_inventory = defaultdict(float)
            for k, v in clone_inventory.items():
                original_inventory[k] = v
            
            clone = {
                'parent': current_planet,
                'offset': offset,
                'aspect_pct': aspect_pct,
                'L': clone_L,
                'inventory': clone_inventory,
                'original_inventory': original_inventory,
                'wasted_inventory': defaultdict(float),
                'debt': clone_debt,
                'type': clone_type
            }
            clones.append(clone)
            all_initial_clones.append({'parent': current_planet, 'offset': offset, 'L': clone_L})
        
        # ---- JUPITER POISON LOGIC (applied before interaction cycle) ----
        if current_planet == 'Jupiter':
            _jp_sign = planet_sign_map.get('Jupiter', '')
            _jp_L = phase5_data['Jupiter']['L']
            _jp_inv = phase5_data['Jupiter']['p5_inventory']
            _jp_current_val = _jp_inv.get('Jupiter', 0.0)

            jupiter_poison_multiplier = 0.0
            jupiter_poison_case = None
            _jp_poison_notes = []  # diagnostic log

            if _jp_current_val <= 0.001:
                _jp_poison_notes.append("[SKIP] Jupiter has no currency to poison (value={:.2f})".format(_jp_current_val))
            else:
                _jp_poison_notes.append("[INFO] Jupiter currency available: {:.2f}".format(_jp_current_val))
                _jp_poison_notes.append("[INFO] Jupiter sign: {}, L: {:.2f}°".format(_jp_sign, _jp_L))

                # --- Helper: check no malefic planet or malefic clone within 22° of Jupiter ---
                def _jp_malefic_free_zone():
                    _malefic_check_planets = ['Saturn', 'Mars', 'Rahu']
                    # Check Ketu if it holds bad currency
                    if phase5_data['Ketu']['p5_inventory'].get('Bad Ketu', 0.0) > 0.001:
                        _malefic_check_planets.append('Ketu')
                        _jp_poison_notes.append("[INFO] Ketu has Bad Ketu currency -> included in malefic check")
                    # Check Moon if malefic
                    if is_moon_malefic_p5():
                        _malefic_check_planets.append('Moon')
                        _jp_poison_notes.append("[INFO] Moon is malefic -> included in malefic check")
                    for _mp in _malefic_check_planets:
                        _mp_L = phase5_data[_mp]['L']
                        _md = abs(_jp_L - _mp_L)
                        if _md > 180: _md = 360 - _md
                        if _md < 22:
                            _jp_poison_notes.append("[FAIL] Malefic-free zone: {} is {:.1f}° from Jupiter (< 22°)".format(_mp, _md))
                            return False
                    # Check malefic virtual clones (from Saturn/Mars already created)
                    for _cl in all_leftover_clones:
                        if _cl['parent'] in ['Saturn', 'Mars']:
                            _cd = abs(_jp_L - _cl['L'])
                            if _cd > 180: _cd = 360 - _cd
                            if _cd < 22:
                                # Check logic: stop only if > 5% bad currency
                                _cl_inv = _cl['inventory']
                                _cl_total_val = sum(v for v in _cl_inv.values() if v > 0.001)
                                _cl_bad_val = sum(v for k, v in _cl_inv.items() if v > 0.001 and 'Bad' in k)
                                _cl_bad_pct = (_cl_bad_val / _cl_total_val * 100.0) if _cl_total_val > 0.001 else 0.0
                                
                                if _cl_bad_pct > 5.0:
                                    _jp_poison_notes.append("[FAIL] Malefic-free zone: Clone({}_H{}) is {:.1f}° from Jupiter (< 22°) with {:.1f}% Bad Currency (>5%)".format(_cl['parent'], _cl['offset'], _cd, _cl_bad_pct))
                                    return False
                                else:
                                    _jp_poison_notes.append("[INFO] Malefic clone ({}_H{}) near ({:.1f}°) but Bad% {:.1f} <= 5% -> Ignored".format(_cl['parent'], _cl['offset'], _cd, _cl_bad_pct))

                    _jp_poison_notes.append("[PASS] Malefic-free zone: No malefic planet or clone (with >5% bad) within 22° of Jupiter")
                    return True

                # --- Case A: Jupiter-Venus Poisoning ---
                _case_a_multiplier = 0.0
                _case_a_signs = {'Sagittarius', 'Pisces', 'Libra', 'Taurus', 'Cancer'}
                _jp_in_parivarthana = 'Jupiter' in parivardhana_map
                _jp_case_a_sign_ok = (_jp_sign in _case_a_signs) or _jp_in_parivarthana
                _jp_poison_notes.append("--- Case A: Jupiter-Venus Poison ---")
                _jp_poison_notes.append("[CHECK] Jupiter sign '{}' in {{Sag,Pis,Lib,Tau,Can}}: {} | Parivarthana: {}".format(_jp_sign, _jp_sign in _case_a_signs, _jp_in_parivarthana))
                _jp_poison_notes.append("[CHECK] Jupiter Case A sign OK (sign or parivarthana): {}".format(_jp_case_a_sign_ok))
                _jp_poison_notes.append("[INFO] Venus sign check: SKIPPED (not required)")
                if _jp_case_a_sign_ok:
                    _venus_L = phase5_data['Venus']['L']
                    _jv_diff = abs(_jp_L - _venus_L)
                    if _jv_diff > 180: _jv_diff = 360 - _jv_diff
                    _jv_gap = int(_jv_diff)
                    _jp_poison_notes.append("[CHECK] Jupiter-Venus distance: {:.1f}° (gap={}) <= 28°: {}".format(_jv_diff, _jv_gap, _jv_gap <= 28))
                    if _jv_gap <= 28:
                        _mfz_a = _jp_malefic_free_zone()
                        if _mfz_a:
                            _cap_pct_a = max(50.0, 100.0 - (_jv_gap * (50.0 / 22.0))) if _jv_gap <= 22 else 50.0
                            _case_a_multiplier = (_cap_pct_a / 100.0) * 0.5
                            _jp_poison_notes.append("[PASS] Case A: gap={}°, cap_pct={:.1f}%, multiplier(cap*50%)={:.2%}".format(_jv_gap, _cap_pct_a, _case_a_multiplier))
                        else:
                            _jp_poison_notes.append("[FAIL] Case A: Malefic-free zone check failed")
                    else:
                        _jp_poison_notes.append("[FAIL] Case A: Jupiter-Venus too far apart ({}° > 28°)".format(_jv_gap))
                else:
                    _jp_poison_notes.append("[FAIL] Case A: Jupiter sign condition not met")

                # --- Case B: Jupiter-Moon Poisoning ---
                _case_b_multiplier = 0.0
                _case_b_signs = {'Sagittarius', 'Pisces', 'Cancer'}
                _jp_case_b_sign_ok = (_jp_sign in _case_b_signs) or _jp_in_parivarthana

                _jp_poison_notes.append("--- Case B: Jupiter-Moon Poison ---")
                _jp_poison_notes.append("[CHECK] Jupiter sign '{}' in {{Sag,Pis,Can}}: {} | Parivarthana: {}".format(_jp_sign, _jp_sign in _case_b_signs, _jp_in_parivarthana))
                _jp_poison_notes.append("[CHECK] Jupiter Case B sign OK (sign or parivarthana): {}".format(_jp_case_b_sign_ok))
                _jp_poison_notes.append("[INFO] Moon sign check: SKIPPED (not required)")

                if _jp_case_b_sign_ok:
                    # Moon phase check: waxing with >50% good OR waning with <2% bad
                    _moon_good_pct = planet_data['Moon'].get('moon_good_pct', 0)
                    _moon_bad_pct = planet_data['Moon'].get('moon_bad_pct', 0)
                    _moon_is_waxing = (paksha == 'Shukla') or (moon_phase_name == 'Purnima')
                    _moon_phase_ok = (_moon_is_waxing and _moon_good_pct > 50) or (not _moon_is_waxing and _moon_bad_pct < 10)
                    _jp_poison_notes.append("[CHECK] Moon phase: paksha={}, waxing={}, good_pct={}, bad_pct={}".format(paksha, _moon_is_waxing, _moon_good_pct, _moon_bad_pct))
                    _jp_poison_notes.append("[CHECK] Moon phase OK (waxing>50% or waning<10% bad): {}".format(_moon_phase_ok))

                    # Moon purity check: Bad Moon < 2% of Moon's total inventory
                    _moon_inv = phase5_data['Moon']['p5_inventory']
                    _moon_bad_currency = _moon_inv.get('Bad Moon', 0.0)
                    _moon_total_inv = sum(abs(v) for v in _moon_inv.values())
                    _moon_pure = _moon_total_inv < 0.001 or (_moon_bad_currency / _moon_total_inv) < 0.02
                    _moon_bad_ratio = (_moon_bad_currency / _moon_total_inv * 100) if _moon_total_inv > 0.001 else 0.0
                    _jp_poison_notes.append("[CHECK] Moon purity: Bad Moon={:.2f}, Total={:.2f}, Bad ratio={:.1f}% < 2%: {}".format(_moon_bad_currency, _moon_total_inv, _moon_bad_ratio, _moon_pure))

                    if _moon_phase_ok and _moon_pure:
                        _moon_L = phase5_data['Moon']['L']
                        _jm_diff = abs(_jp_L - _moon_L)
                        if _jm_diff > 180: _jm_diff = 360 - _jm_diff
                        _jm_gap = int(_jm_diff)
                        _jp_poison_notes.append("[CHECK] Jupiter-Moon distance: {:.1f}° (gap={}) <= 28°: {}".format(_jm_diff, _jm_gap, _jm_gap <= 28))
                        if _jm_gap <= 28:
                            _mfz_b = _jp_malefic_free_zone()
                            if _mfz_b:
                                _cap_pct_b = max(50.0, 100.0 - (_jm_gap * (50.0 / 22.0))) if _jm_gap <= 22 else 50.0
                                _case_b_multiplier = _cap_pct_b / 100.0
                                _jp_poison_notes.append("[PASS] Case B: gap={}°, cap_pct={:.1f}%, multiplier={:.2%}".format(_jm_gap, _cap_pct_b, _case_b_multiplier))
                            else:
                                _jp_poison_notes.append("[FAIL] Case B: Malefic-free zone check failed")
                        else:
                            _jp_poison_notes.append("[FAIL] Case B: Jupiter-Moon too far apart ({}° > 28°)".format(_jm_gap))
                    else:
                        if not _moon_phase_ok:
                            _jp_poison_notes.append("[FAIL] Case B: Moon phase condition not met")
                        if not _moon_pure:
                            _jp_poison_notes.append("[FAIL] Case B: Moon purity condition not met (Bad Moon ratio {:.1f}% >= 2%)".format(_moon_bad_ratio))
                else:
                    _jp_poison_notes.append("[FAIL] Case B: Jupiter sign condition not met")

                # --- Pick the highest multiplier if both qualify ---
                _jp_poison_notes.append("--- Result ---")
                _jp_poison_notes.append("[INFO] Case A multiplier: {:.2%}, Case B multiplier: {:.2%}".format(_case_a_multiplier, _case_b_multiplier))
                if _case_a_multiplier > 0.001 or _case_b_multiplier > 0.001:
                    if _case_a_multiplier >= _case_b_multiplier:
                        jupiter_poison_multiplier = _case_a_multiplier
                        jupiter_poison_case = 'CaseA_Venus'
                    else:
                        jupiter_poison_multiplier = _case_b_multiplier
                        jupiter_poison_case = 'CaseB_Moon'
                    _jp_poison_notes.append("[APPLIED] Winner: {} with multiplier {:.2%}".format(jupiter_poison_case, jupiter_poison_multiplier))
                else:
                    _jp_poison_notes.append("[NOT APPLIED] Neither case qualified — no poison applied")

                # --- Apply poison to Jupiter's own inventory ---
                if jupiter_poison_multiplier > 0.001:
                    _poison_amount = jupiter_poison_multiplier * _jp_current_val
                    _jp_inv['Jupiter'] = _jp_current_val - _poison_amount
                    _jp_inv['Jupiter Poison'] = _jp_inv.get('Jupiter Poison', 0.0) + _poison_amount
                    _jp_poison_notes.append("[ACTION] Jupiter own: {:.2f} -> Jupiter[{:.2f}] + Poison[{:.2f}]".format(_jp_current_val, _jp_current_val - _poison_amount, _poison_amount))

                    # --- Apply poison to all Jupiter clones ---
                    for _jcl in clones:
                        _jcl_jup_val = _jcl['inventory'].get('Jupiter', 0.0)
                        if _jcl_jup_val > 0.001:
                            _cl_poison = jupiter_poison_multiplier * _jcl_jup_val
                            _jcl['inventory']['Jupiter'] = _jcl_jup_val - _cl_poison
                            _jcl['inventory']['Jupiter Poison'] = _jcl['inventory'].get('Jupiter Poison', 0.0) + _cl_poison
                            # Also update original_inventory so Steps 2 & 3 pull poisoned values
                            _jcl_orig_val = _jcl['original_inventory'].get('Jupiter', 0.0)
                            if _jcl_orig_val > 0.001:
                                _cl_orig_poison = jupiter_poison_multiplier * _jcl_orig_val
                                _jcl['original_inventory']['Jupiter'] = _jcl_orig_val - _cl_orig_poison
                                _jcl['original_inventory']['Jupiter Poison'] = _jcl['original_inventory'].get('Jupiter Poison', 0.0) + _cl_orig_poison
                            _jp_poison_notes.append("[ACTION] Clone(H{}): {:.2f} -> Jupiter[{:.2f}] + Poison[{:.2f}]".format(_jcl['offset'], _jcl_jup_val, _jcl_jup_val - _cl_poison, _cl_poison))

        # Part B: The Interaction Cycle
        # MODIFICATION 2: Reordered - Step 1 is Active Pulling, Step 2 is Real Malefics Pull
        p5_cycle_limit = 500
        p5_cycles = 0
        
        while p5_cycles < p5_cycle_limit:
            p5_cycles += 1
            p5_something_happened = False
            
            for clone in clones:
                clone_L = clone['L']
                
                # MODIFICATION 2: NEW Step 1 - Virtual Malefic Clones Pull (Active Pulling)
                if clone['type'] == 'Active' and clone['debt'] < -0.001:
                    clone_parent = clone['parent']
                    
                    targets_by_proximity = []
                    for target_planet in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
                        target_L = phase5_data[target_planet]['L']
                        diff = abs(clone_L - target_L)
                        if diff > 180: diff = 360 - diff
                        gap = int(diff)
                        if gap <= 22:
                            targets_by_proximity.append({
                                'planet': target_planet,
                                'gap': gap,
                                'L': target_L
                            })
                    
                    targets_by_proximity.sort(key=lambda x: x['gap'])
                    
                    for target_info in targets_by_proximity:
                        if clone['debt'] >= -0.001:
                            break
                        
                        target_planet = target_info['planet']
                        gap = target_info['gap']
                        
                        target_vol = phase5_data[target_planet]['volume']
                        cap_pct = mix_dict.get(gap, 0)
                        max_allowed_pull = target_vol * (cap_pct / 100.0)
                        
                        tracker_key = f"clone_{clone['parent']}_{clone['offset']}_pulled_from_{target_planet}"
                        already_pulled = clone.get(tracker_key, 0.0)
                        remaining_capacity = max_allowed_pull - already_pulled
                        
                        if remaining_capacity <= 0.001:
                            continue
                        
                        banned_currency = f"Good {clone_parent}"
                        
                        target_inv = phase5_data[target_planet]['p5_inventory']
                        
                        good_currencies = []
                        bad_currencies = []
                        for k, v in target_inv.items():
                            if v > 0.001 and k != banned_currency:
                                score = get_p5_currency_rank_score(k)
                                if is_good_currency(k):
                                    good_currencies.append({'key': k, 'value': v, 'score': score})
                                else:
                                    bad_currencies.append({'key': k, 'value': v, 'score': score})
                        
                        good_currencies.sort(key=lambda x: -x['score'])
                        bad_currencies.sort(key=lambda x: -x['score'])
                        
                        for curr in good_currencies:
                            if clone['debt'] >= -0.001:
                                break
                            if remaining_capacity <= 0.001:
                                break
                            
                            avail = target_inv[curr['key']]
                            if avail <= 0.001:
                                continue
                            
                            needed_debt = abs(clone['debt'])
                            take = min(1.0, needed_debt, avail, remaining_capacity)
                            
                            if take > 0.001:
                                phase5_data[target_planet]['p5_inventory'][curr['key']] -= take
                                phase5_data[target_planet]['p5_current_debt'] -= take
                                
                                # MODIFICATION 2: Add to clone's inventory AND mark as wasted
                                clone['inventory'][curr['key']] = clone['inventory'].get(curr['key'], 0.0) + take
                                clone['wasted_inventory'][curr['key']] = clone['wasted_inventory'].get(curr['key'], 0.0) + take
                                
                                clone['debt'] += take
                                
                                clone[tracker_key] = already_pulled + take
                                remaining_capacity -= take
                                already_pulled += take
                                p5_something_happened = True
                        
                        good_still_available = any(
                            target_inv.get(c['key'], 0) > 0.001 
                            for c in good_currencies 
                            if c['key'] != banned_currency
                        )
                        
                        if clone['debt'] < -0.001 and not good_still_available:
                            for curr in bad_currencies:
                                if clone['debt'] >= -0.001:
                                    break
                                if remaining_capacity <= 0.001:
                                    break
                                
                                avail = target_inv[curr['key']]
                                if avail <= 0.001:
                                    continue
                                
                                needed_debt = abs(clone['debt'])
                                take = min(1.0, needed_debt, avail, remaining_capacity)
                                
                                if take > 0.001:
                                    phase5_data[target_planet]['p5_inventory'][curr['key']] -= take
                                    phase5_data[target_planet]['p5_current_debt'] -= take
                                    if 'Bad' in curr['key']:
                                        phase5_data[target_planet]['bad_inv'] -= take
                                    
                                    clone['inventory'][curr['key']] = clone['inventory'].get(curr['key'], 0.0) + take
                                    clone['wasted_inventory'][curr['key']] = clone['wasted_inventory'].get(curr['key'], 0.0) + take
                                    
                                    clone['debt'] += take
                                    
                                    clone[tracker_key] = already_pulled + take
                                    remaining_capacity -= take
                                    already_pulled += take
                                    p5_something_happened = True
                
                # MODIFICATION 2: NEW Step 2 - Real Malefics Pull (from Clone's original inventory only)
                total_original_remaining = 0.0
                for k in clone['original_inventory'].keys():
                    taken_key = f'taken_from_original_{k}'
                    taken = clone.get(taken_key, 0.0)
                    remaining = clone['original_inventory'][k] - taken
                    if remaining > 0.001:
                        total_original_remaining += remaining
                
                if total_original_remaining > 0.001:
                    real_malefics = list(P5_MALEFIC_DEBTOR_RANK)
                    if is_moon_malefic_p5() and 'Moon' not in real_malefics:
                        moon_vol = phase5_data['Moon']['volume']
                        if moon_vol > 0:
                            bad_pct = (phase5_data['Moon']['bad_inv'] / moon_vol) * 100
                            if bad_pct > 25:
                                mars_idx = real_malefics.index('Mars')
                                real_malefics.insert(mars_idx, 'Moon')
                            else:
                                mars_idx = real_malefics.index('Mars')
                                real_malefics.insert(mars_idx + 1, 'Moon')
                    
                    for malefic in real_malefics:
                        if phase5_data[malefic]['p5_current_debt'] >= -0.001:
                            continue
                        
                        malefic_L = phase5_data[malefic]['L']
                        diff = abs(malefic_L - clone_L)
                        if diff > 180: diff = 360 - diff
                        gap = int(diff)
                        if gap > 22:
                            continue
                        
                        total_original_vol = sum(clone['original_inventory'].values())
                        cap_pct = mix_dict.get(gap, 0)
                        max_allowed_pull = total_original_vol * (cap_pct / 100.0)
                        
                        tracker_key = f"p5_pulled_from_clone_{clone['parent']}_{clone['offset']}"
                        already_pulled = phase5_data[malefic].get(tracker_key, 0.0)
                        remaining_capacity = max_allowed_pull - already_pulled
                        
                        if remaining_capacity <= 0.001:
                            continue
                        
                        available_currencies = []
                        for k, orig_v in clone['original_inventory'].items():
                            taken_key = f'taken_from_original_{k}'
                            taken = clone.get(taken_key, 0.0)
                            remaining_original = orig_v - taken
                            if remaining_original > 0.001 and is_good_currency(k):
                                score = get_p5_currency_rank_score(k)
                                available_currencies.append({'key': k, 'remaining': remaining_original, 'score': score})
                        
                        available_currencies.sort(key=lambda x: -x['score'])
                        
                        for curr in available_currencies:
                            if phase5_data[malefic]['p5_current_debt'] >= -0.001:
                                break
                            if remaining_capacity <= 0.001:
                                break
                            
                            taken_key = f"taken_from_original_{curr['key']}"
                            taken = clone.get(taken_key, 0.0)
                            remaining_original = clone['original_inventory'][curr['key']] - taken
                            
                            if remaining_original <= 0.001:
                                continue
                            
                            needed_debt = abs(phase5_data[malefic]['p5_current_debt'])
                            take = min(1.0, needed_debt, remaining_original, remaining_capacity)
                            
                            if take > 0.001:
                                clone[taken_key] = taken + take
                                clone['inventory'][curr['key']] -= take
                                
                                phase5_data[malefic]['p5_inventory'][curr['key']] += take
                                phase5_data[malefic]['p5_current_debt'] += take
                                
                                phase5_data[malefic][tracker_key] = already_pulled + take
                                remaining_capacity -= take
                                already_pulled += take
                                p5_something_happened = True
                
                # Step 3: Real Benefics Pull (from Clone's original inventory only)
                total_original_remaining = 0.0
                for k in clone['original_inventory'].keys():
                    taken_key = f'taken_from_original_{k}'
                    taken = clone.get(taken_key, 0.0)
                    remaining = clone['original_inventory'][k] - taken
                    if remaining > 0.001:
                        total_original_remaining += remaining
                
                if total_original_remaining > 0.001:
                    real_benefics = list(P5_STANDARD_BENEFICS)
                    if not is_moon_malefic_p5():
                        real_benefics.append('Moon')
                    
                    def get_benefic_debt_pct_p5(p):
                        vol = phase5_data[p]['volume']
                        debt = abs(phase5_data[p]['p5_current_debt'])
                        if vol > 0:
                            return (debt / vol) * 100
                        return 0
                    
                    real_benefics_sorted = sorted(real_benefics, key=lambda p: -get_benefic_debt_pct_p5(p))
                    
                    for benefic in real_benefics_sorted:
                        if phase5_data[benefic]['p5_current_debt'] >= -0.001:
                            continue
                        
                        benefic_L = phase5_data[benefic]['L']
                        diff = abs(benefic_L - clone_L)
                        if diff > 180: diff = 360 - diff
                        gap = int(diff)
                        if gap > 22:
                            continue
                        
                        total_original_vol = sum(clone['original_inventory'].values())
                        cap_pct = mix_dict.get(gap, 0)
                        max_allowed_pull = total_original_vol * (cap_pct / 100.0)
                        
                        tracker_key = f"p5_benefic_pulled_from_clone_{clone['parent']}_{clone['offset']}"
                        already_pulled = phase5_data[benefic].get(tracker_key, 0.0)
                        remaining_capacity = max_allowed_pull - already_pulled
                        
                        if remaining_capacity <= 0.001:
                            continue
                        
                        available_currencies = []
                        for k, orig_v in clone['original_inventory'].items():
                            taken_key = f'taken_from_original_{k}'
                            taken = clone.get(taken_key, 0.0)
                            remaining_original = orig_v - taken
                            if remaining_original > 0.001 and is_good_currency(k):
                                score = get_p5_currency_rank_score(k)
                                available_currencies.append({'key': k, 'remaining': remaining_original, 'score': score})
                        
                        available_currencies.sort(key=lambda x: -x['score'])
                        
                        for curr in available_currencies:
                            if phase5_data[benefic]['p5_current_debt'] >= -0.001:
                                break
                            if remaining_capacity <= 0.001:
                                break
                            
                            taken_key = f"taken_from_original_{curr['key']}"
                            taken = clone.get(taken_key, 0.0)
                            remaining_original = clone['original_inventory'][curr['key']] - taken
                            
                            if remaining_original <= 0.001:
                                continue
                            
                            needed_debt = abs(phase5_data[benefic]['p5_current_debt'])
                            take = min(1.0, needed_debt, remaining_original, remaining_capacity)
                            
                            if take > 0.001:
                                clone[taken_key] = taken + take
                                clone['inventory'][curr['key']] -= take
                                
                                phase5_data[benefic]['p5_inventory'][curr['key']] += take
                                phase5_data[benefic]['p5_current_debt'] += take
                                
                                phase5_data[benefic][tracker_key] = already_pulled + take
                                remaining_capacity -= take
                                already_pulled += take
                                p5_something_happened = True
            
            if not p5_something_happened:
                break
        
        # Part C: Logging & Disposal
        for clone in clones:
            inv_parts = []
            for k, v in clone['inventory'].items():
                if v > 0.001:
                    inv_parts.append(f"{k}[{v:.2f}]")
            inv_str = ", ".join(inv_parts) if inv_parts else "-"
            
            if abs(clone['debt']) < 0.01:
                debt_str = "0.00"
            else:
                debt_str = f"{clone['debt']:.2f}"
            
            leftover_aspects.append([
                clone['parent'],
                clone['offset'],
                f"{clone['L']:.2f}",
                inv_str,
                debt_str
            ])
            all_leftover_clones.append(clone)
    
    # ---- JUPITER POISON POST-SHARING DEBT APPLICATION ----
    # After all sharing is done, add -2 debt per unit of Jupiter Poison held
    # This converts Jupiter Poison from appearing good during sharing to being penalised
    for _jp_p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        _jp_poison_held = phase5_data[_jp_p]['p5_inventory'].get('Jupiter Poison', 0.0)
        if _jp_poison_held > 0.001:
            phase5_data[_jp_p]['p5_current_debt'] -= 2.0 * _jp_poison_held

    # Apply debt to leftover clones holding Jupiter Poison
    for _jp_cl in all_leftover_clones:
        _jp_cl_poison = _jp_cl['inventory'].get('Jupiter Poison', 0.0)
        if _jp_cl_poison > 0.001:
            _jp_cl['debt'] -= 2.0 * _jp_cl_poison

    # ---- KETU ALONE & UNASPECTED CHECK ----
    # If Ketu is not conjuncted or aspected by any planet (within 22 degrees)
    # and resides alone in Gemini, Leo, Scorpio, or Aquarius => add -25 Bad Ketu and -25 debt
    _ketu_lonely_signs = {'Gemini', 'Leo', 'Scorpio', 'Aquarius'}
    _ketu_sign = planet_sign_map.get('Ketu', '')
    if _ketu_sign in _ketu_lonely_signs:
        _ketu_L = phase5_data['Ketu']['L']
        _ketu_is_alone = True
        _all_planets_for_ketu_check = ['Sun', 'Moon', 'Mars', 'Mercury', 'Jupiter', 'Venus', 'Saturn', 'Rahu']
        for _chk_p in _all_planets_for_ketu_check:
            _chk_L = phase5_data[_chk_p]['L']
            _raw_diff = abs(_ketu_L - _chk_L)
            if _raw_diff > 180:
                _raw_diff = 360 - _raw_diff
            if _raw_diff < 22:
                _ketu_is_alone = False
                break
        # Also check aspect clones landing near Ketu (within 22 degrees)
        if _ketu_is_alone:
            for _cl in all_leftover_clones:
                _cl_L = _cl['L']
                _raw_diff = abs(_ketu_L - _cl_L)
                if _raw_diff > 180:
                    _raw_diff = 360 - _raw_diff
                if _raw_diff < 22:
                    _ketu_is_alone = False
                    break
        if _ketu_is_alone:
            phase5_data['Ketu']['p5_inventory']['Bad Ketu'] = phase5_data['Ketu']['p5_inventory'].get('Bad Ketu', 0.0) + 25.0
            phase5_data['Ketu']['p5_current_debt'] -= 25.0

    # Format Phase 5 Output
    phase5_rows = []
    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        inv = phase5_data[p]['p5_inventory']
        parts = []
        own_keys = [p, f"Good {p}", f"Bad {p}"]
        if p == 'Moon': own_keys = ["Good Moon", "Bad Moon"]
        for k in own_keys:
            if k in inv and inv[k] > 0.001: parts.append(f"{k}[{inv[k]:.2f}]")
        for k, v in inv.items():
            if k not in own_keys and v > 0.001: parts.append(f"{k}[{v:.2f}]")
        phase5_data[p]['currency_p5'] = ", ".join(parts) if parts else "-"
        d_val = phase5_data[p]['p5_current_debt']
        if abs(d_val) < 0.01: phase5_data[p]['debt_p5'] = "0.00"
        else: phase5_data[p]['debt_p5'] = f"{d_val:.2f}"
    
    # Jupiter Poison penalty multipliers
    poisonpenality = 3    # used in HPS (aspect + occupant) and NPS
    poisonpenality_1 = 2  # used in Phase 5 Net Currency Score

    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        d_p5 = phase5_data[p]
        inv = phase5_data[p]['p5_inventory']
        # Jupiter Poison treated as bad for net score after sharing (poisonpenality_1 = 2x)
        net_score = sum(
            (-poisonpenality_1*v if k == 'Jupiter Poison' else (v if is_good_currency(k) else -v))
            for k, v in inv.items()
        )
        phase5_rows.append([p, d_p5['currency_p5'], d_p5['debt_p5'], f"{net_score:.2f}"])
    
    df_phase5 = pd.DataFrame(phase5_rows, columns=['Planet', 'Currency [Phase 5]', 'Debt [Phase 5]', 'Net Currency Score'])
    
    df_leftover_aspects = pd.DataFrame(leftover_aspects, columns=['Source Planet', 'Aspect Angle', 'Position', 'Remaining Inventory', 'Final Debt'])
    
    # JUPITER POISON DIAGNOSTIC NOTES
    _jp_diag_rows = []
    if not _jp_poison_notes:
        _jp_poison_notes = ["[INFO] No Jupiter Poison conditions were evaluated"]
    for _note_idx, _note_line in enumerate(_jp_poison_notes, 1):
        _jp_diag_rows.append([_note_idx, _note_line])
    # Post-sharing poison debt summary
    _jp_debt_notes = []
    for _jdp in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        _jdp_poison = phase5_data[_jdp]['p5_inventory'].get('Jupiter Poison', 0.0)
        if _jdp_poison > 0.001:
            _jp_debt_notes.append("{}: holds {:.2f} Jupiter Poison -> debt penalty: {:.2f}".format(_jdp, _jdp_poison, -2.0 * _jdp_poison))
    for _jcl_d in all_leftover_clones:
        _jcl_d_poison = _jcl_d['inventory'].get('Jupiter Poison', 0.0)
        if _jcl_d_poison > 0.001:
            _jp_debt_notes.append("Clone({}_H{}): holds {:.2f} Jupiter Poison -> debt penalty: {:.2f}".format(_jcl_d['parent'], _jcl_d['offset'], _jcl_d_poison, -2.0 * _jcl_d_poison))
    if _jp_debt_notes:
        _jp_diag_rows.append([len(_jp_diag_rows) + 1, "--- Post-Sharing Debt Application ---"])
        for _dbn in _jp_debt_notes:
            _jp_diag_rows.append([len(_jp_diag_rows) + 1, _dbn])
    else:
        _jp_diag_rows.append([len(_jp_diag_rows) + 1, "[INFO] No planet or clone holds Jupiter Poison after sharing"])
    df_jupiter_poison_notes = pd.DataFrame(_jp_diag_rows, columns=['#', 'Jupiter Poison Diagnostic'])

    # CREATE HOUSE RESERVES DATAFRAME
    reserve_rows = []
    for sign_name in sign_names:
        currency_dict = house_reserves[sign_name]
        
        currency_parts = []
        for currency_type, amount in currency_dict.items():
            if amount > 0.001:
                currency_parts.append(f"{currency_type}[{amount:.2f}]")
        
        reserve_str = ", ".join(currency_parts) if currency_parts else "-"
        reserve_rows.append([sign_name, reserve_str])
    
    df_house_reserves = pd.DataFrame(reserve_rows, columns=['House Sign', 'Unutilized Bonus Points'])

    # ---- HOUSE POINTS ANALYSIS (Pre-calc for KHS) ----
    aspect_score   = {s: 0.0 for s in sign_names}
    aspect_sources = {s: [] for s in sign_names}
    occupant_score = {s: 0.0 for s in sign_names}
    occupant_notes = {s: [] for s in sign_names}

    # ---- NATURAL PLANETARY RELATIONSHIPS (Moved here for use in Aspect Logic) ----
    NATURAL_FRIENDSHIPS = {
        'Sun': {'Friends': ['Moon', 'Mars', 'Jupiter'], 'Neutral': ['Mercury'], 'Enemies': ['Venus', 'Saturn']},
        'Moon': {'Friends': ['Sun', 'Mercury'], 'Neutral': ['Mars', 'Jupiter', 'Venus', 'Saturn'], 'Enemies': []},
        # MODIFIED: Venus & Saturn moved to Enemies
        'Mars': {'Friends': ['Sun', 'Moon', 'Jupiter'], 'Neutral': [], 'Enemies': ['Venus', 'Saturn', 'Mercury']}, 
        'Mercury': {'Friends': ['Sun', 'Venus'], 'Neutral': ['Mars', 'Jupiter', 'Saturn'], 'Enemies': ['Moon']},
        'Jupiter': {'Friends': ['Sun', 'Moon', 'Mars'], 'Neutral': ['Saturn'], 'Enemies': ['Mercury', 'Venus']},
        'Venus': {'Friends': ['Mercury', 'Saturn'], 'Neutral': ['Mars', 'Jupiter'], 'Enemies': ['Sun', 'Moon']},
        # MODIFIED: Jupiter moved to Enemies
        'Saturn': {'Friends': ['Mercury', 'Venus'], 'Neutral': [], 'Enemies': ['Sun', 'Moon', 'Mars', 'Jupiter']},
        'Rahu': {'Friends': ['Mercury', 'Venus', 'Saturn'], 'Neutral': ['Jupiter'], 'Enemies': ['Sun', 'Moon', 'Mars']},
        'Ketu': {'Friends': ['Mars', 'Venus', 'Saturn'], 'Neutral': ['Mercury', 'Jupiter'], 'Enemies': ['Sun', 'Moon']}
    }

    def check_friendship(planet, target):
        """Check relationship of planet toward target using NATURAL_FRIENDSHIPS."""
        if planet == target:
            return 'Friend'
        rels = NATURAL_FRIENDSHIPS.get(planet, {})
        if target in rels.get('Friends', []):
            return 'Friend'
        elif target in rels.get('Enemies', []):
            return 'Enemy'
        else:
            return 'Neutral'

    hp_static_malefics = {'Saturn', 'Mars', 'Sun', 'Rahu'}
    hp_static_benefics = {'Jupiter', 'Venus', 'Mercury'}

    def _hp_is_malefic(planet_name):
        if planet_name in hp_static_malefics:
            return True
        if planet_name == 'Moon':
            return phase5_data['Moon']['bad_inv'] > 0.001
        if planet_name == 'Ketu':
            return phase5_data['Ketu']['p5_inventory'].get('Bad Ketu', 0.0) > 0.001
        return False

    _own_good_key = {
        'Saturn': 'Good Saturn', 'Mars': 'Good Mars', 'Sun': 'Good Sun',
        'Rahu': 'Good Rahu', 'Ketu': 'Good Ketu', 'Moon': 'Good Moon'
    }

    lagna_sign_hp = get_sign(lagna_sid)
    lagna_lord = get_sign_lord(lagna_sign_hp)
    is_malefic_lagna_lord = (
        lagna_lord in ('Sun', 'Mars', 'Saturn')
        or (lagna_lord == 'Moon' and phase5_data['Moon']['bad_inv'] > 0.001)
    )

    for clone in all_leftover_clones:
        parent = clone['parent']
        parent_L = phase5_data[parent]['L']
        target_lon = (parent_L + (clone['offset'] - 1) * 30) % 360
        target_sign = get_sign(target_lon)

        if _hp_is_malefic(parent):
            if clone['debt'] < -0.001:
                if parent == lagna_lord and is_malefic_lagna_lord and target_sign == lagna_sign_hp:
                    penalty = abs(clone['debt']) / 2.0
                    aspect_score[target_sign] -= penalty
                    aspect_sources[target_sign].append(f"{parent}(Lagna Lord Debt/2 [Lagna])")
                else:
                    penalty = abs(clone['debt'])
                    aspect_score[target_sign] -= penalty
                    aspect_sources[target_sign].append(f"{parent}(Malefic Debt)")
            own_key = _own_good_key.get(parent)
            if own_key:
                own_val = clone['inventory'].get(own_key, 0.0)
                if own_val > 0.001:
                    # NEW LOGIC: Mars/Saturn check specifically for Enemy House
                    is_enemy_house = False
                    if parent in ['Mars', 'Saturn']:
                         # Special rule: Mars in Leo is never negative for enemy houses
                         if parent == 'Mars' and planet_data['Mars']['sign'] == 'Leo':
                             is_enemy_house = False
                         else:
                             target_lord = get_sign_lord(target_sign)
                             relation = check_friendship(parent, target_lord)
                             if relation == 'Enemy':
                                 is_enemy_house = True
                    
                    if is_enemy_house:
                        aspect_score[target_sign] -= own_val
                        aspect_sources[target_sign].append(f"{parent}(Own Good converted to Neg [Enemy House])")
                    else:
                        bonus = own_val
                        aspect_score[target_sign] += bonus
                        aspect_sources[target_sign].append(f"{parent}(Own Good)")
        else:
            good_total = 0.0
            for c_key, c_val in clone['inventory'].items():
                if c_val > 0.001 and is_good_currency(c_key):
                    good_total += c_val
            # Jupiter Poison: treat as bad in HPS (subtract from benefic bonus)
            _hp_cl_poison = clone['inventory'].get('Jupiter Poison', 0.0)
            if _hp_cl_poison > 0.001:
                good_total -= _hp_cl_poison  # remove poison from good
            if good_total > 0.001:
                aspect_score[target_sign] += good_total
                aspect_sources[target_sign].append(f"{parent}(Benefic Bonus)")
            # Subtract Jupiter Poison as a penalty (poisonpenality - 1 extra beyond removal)
            if _hp_cl_poison > 0.001:
                aspect_score[target_sign] -= (poisonpenality - 1) * _hp_cl_poison
                aspect_sources[target_sign].append(f"{parent}(Jupiter Poison Penalty x{poisonpenality})")

    sign_occupants = defaultdict(list)
    for p_name in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        occ_sign = planet_sign_map[p_name]
        sign_occupants[occ_sign].append(p_name)

    # ---- KENDRAADHIBATHYA DOSHA CHECK ----
    # Config: (lagna_sign, planet, planet_sign, penalty)
    _kad_configs = [
        ('Sagittarius', 'Jupiter', 'Pisces', -150.0),
        ('Gemini', 'Mercury', 'Virgo', -70.0),
    ]
    _kad_results = {}  # planet -> target_sign (for skipping in occupant loop & gift pot)

    for _kad_lagna, _kad_planet, _kad_target, _kad_penalty in _kad_configs:
        if lagna_sign_hp == _kad_lagna and planet_sign_map.get(_kad_planet) == _kad_target:
            _kad_pL = phase5_data[_kad_planet]['L']
            _kad_mf = True
            _kad_chk = ['Saturn', 'Mars', 'Rahu']
            # Include Ketu only if it holds Bad Ketu currency
            if phase5_data['Ketu']['p5_inventory'].get('Bad Ketu', 0.0) > 0.001:
                _kad_chk.append('Ketu')
            # Include Moon if malefic
            if is_moon_malefic_p5():
                _kad_chk.append('Moon')
            for _km in _kad_chk:
                if _km == _kad_planet:
                    continue  # don't check planet against itself
                _km_L = phase5_data[_km]['L']
                _km_d = abs(_kad_pL - _km_L)
                if _km_d > 180: _km_d = 360 - _km_d
                if _km_d < 22:
                    _kad_mf = False
                    break
            # Also check malefic virtual clones within 22°
            if _kad_mf:
                for _kcl in all_leftover_clones:
                    if _kcl['parent'] in ['Saturn', 'Mars', 'Rahu', 'Ketu']:
                        _kcl_d = abs(_kad_pL - _kcl['L'])
                        if _kcl_d > 180: _kcl_d = 360 - _kcl_d
                        if _kcl_d < 22:
                            _kcl_inv = _kcl['inventory']
                            _kcl_bad = sum(v for k, v in _kcl_inv.items() if v > 0.001 and 'Bad' in k)
                            if _kcl_bad > 0.001:
                                _kad_mf = False
                                break
            if _kad_mf:
                _kad_results[_kad_planet] = _kad_target
                occupant_score[_kad_target] += _kad_penalty
                occupant_notes[_kad_target].append(f"{_kad_planet}({_kad_penalty:.0f} Kendraadhibathya Dosha)")

    # ---- KONA DOSHA CHECK ----
    # Pisces lagna + Jupiter in Cancer: if Jupiter is malefic-free AND Pisces (1st house) is
    # malefic-free (no malefic planet sitting in Pisces, no malefic clone aspecting Pisces),
    # penalize 1st house occupant score by -100. Benefics in Pisces still contribute normally.
    _kona_dosha_active = False
    if lagna_sign_hp == 'Pisces' and planet_sign_map.get('Jupiter') == 'Cancer':
        _kona_jp_L = phase5_data['Jupiter']['L']
        _kona_pass = True

        # Step 1: Check Jupiter is malefic-free (no malefic planet within 22° of Jupiter)
        _kona_mal_list = ['Saturn', 'Mars', 'Rahu']
        if phase5_data['Ketu']['p5_inventory'].get('Bad Ketu', 0.0) > 0.001:
            _kona_mal_list.append('Ketu')
        if is_moon_malefic_p5():
            _kona_mal_list.append('Moon')

        for _km in _kona_mal_list:
            _km_L = phase5_data[_km]['L']
            _km_d = abs(_kona_jp_L - _km_L)
            if _km_d > 180: _km_d = 360 - _km_d
            if _km_d < 22:
                _kona_pass = False
                break

        # Check malefic virtual clones within 22° of Jupiter
        if _kona_pass:
            for _kcl in all_leftover_clones:
                if _kcl['parent'] in ['Saturn', 'Mars', 'Rahu', 'Ketu']:
                    _kcl_d = abs(_kona_jp_L - _kcl['L'])
                    if _kcl_d > 180: _kcl_d = 360 - _kcl_d
                    if _kcl_d < 22:
                        _kcl_bad = sum(v for k, v in _kcl['inventory'].items() if v > 0.001 and 'Bad' in k)
                        if _kcl_bad > 0.001:
                            _kona_pass = False
                            break

        # Step 2: Check Pisces (1st house) is malefic-free
        # 2a: No malefic planet sitting in Pisces
        if _kona_pass:
            for _km in _kona_mal_list:
                if planet_sign_map.get(_km) == 'Pisces':
                    _kona_pass = False
                    break

        # 2b: No malefic clone aspecting (targeting) Pisces
        if _kona_pass:
            for _kcl in all_leftover_clones:
                if _kcl['parent'] in ['Saturn', 'Mars', 'Rahu', 'Ketu']:
                    _kcl_parent_L = phase5_data[_kcl['parent']]['L']
                    _kcl_target_lon = (_kcl_parent_L + (_kcl['offset'] - 1) * 30) % 360
                    _kcl_target_sign = get_sign(_kcl_target_lon)
                    if _kcl_target_sign == 'Pisces':
                        _kcl_bad = sum(v for k, v in _kcl['inventory'].items() if v > 0.001 and 'Bad' in k)
                        if _kcl_bad > 0.001:
                            _kona_pass = False
                            break

        if _kona_pass:
            _kona_dosha_active = True
            occupant_score['Pisces'] += -100.0
            occupant_notes['Pisces'].append("Jupiter(-100 Kona Dosha)")

    for s in sign_names:
        for occ in sign_occupants.get(s, []):
            if occ == 'Rahu':
                continue  # Rahu uses Rahu Score directly, applied after NPS calculation
            # Skip planet's occupant contribution if Kendraadhibathya Dosha is active for it
            if occ in _kad_results and s == _kad_results[occ]:
                occupant_notes[s].append(f"{occ}(Good skipped - Kendraadhibathya Dosha)")
                continue
            inv = phase5_data[occ]['p5_inventory']
            if _hp_is_malefic(occ):
                total_good = sum(v for k, v in inv.items() if v > 0.001 and is_good_currency(k))
                total_bad  = sum(v for k, v in inv.items() if v > 0.001 and 'Bad' in k)
                # Jupiter Poison: treat as bad for HPS occupant scoring (poisonpenality = 3x)
                _hp_occ_poison = inv.get('Jupiter Poison', 0.0)
                if _hp_occ_poison > 0.001:
                    total_good -= _hp_occ_poison
                    total_bad += (poisonpenality - 1) * _hp_occ_poison
                if occ == lagna_lord and is_malefic_lagna_lord and s == lagna_sign_hp:
                    total_bad = total_bad / 2.0
                    net = total_good - total_bad
                    occupant_score[s] += net
                    occupant_notes[s].append(f"{occ}(Good-Bad/2 [LL in Lagna])")
                else:
                    net = total_good - total_bad
                    occupant_score[s] += net
                    occupant_notes[s].append(f"{occ}(Good-Bad)")
            else:
                total_good = sum(v for k, v in inv.items() if v > 0.001 and is_good_currency(k))
                # Jupiter Poison: treat as bad for HPS benefic occupant scoring (poisonpenality = 3x)
                _hp_occ_poison_b = inv.get('Jupiter Poison', 0.0)
                if _hp_occ_poison_b > 0.001:
                    total_good -= _hp_occ_poison_b
                if total_good > 0.001:
                    occupant_score[s] += total_good
                    occupant_notes[s].append(f"{occ}(Benefic Sum)")
                if _hp_occ_poison_b > 0.001:
                    occupant_score[s] -= (poisonpenality - 1) * _hp_occ_poison_b
                    occupant_notes[s].append(f"{occ}(Jupiter Poison Penalty x{poisonpenality})")

    hp_gift_pot_config = {
        'Sagittarius': ('Jupiter', 100),
        'Pisces': ('Jupiter', 80),
        'Libra': ('Venus', 80),
        'Taurus': ('Venus', 60)
    }
    for gift_sign, (gifter, multiplier) in hp_gift_pot_config.items():
        # Skip unused bonus if Kendraadhibathya Dosha is active for this sign
        if gifter in _kad_results and gift_sign == _kad_results[gifter]:
            occupant_notes[gift_sign].append("Unused Bonus skipped (Kendraadhibathya Dosha)")
            continue
        gifter_sthana = planet_data[gifter]['sthana']
        original_pot = multiplier * (gifter_sthana / 100.0)
        unused_reserve = sum(v for v in house_reserves[gift_sign].values() if v > 0.001)
        used_amount = original_pot - unused_reserve
        if used_amount > 0.001:
            occupant_score[gift_sign] -= used_amount
            occupant_notes[gift_sign].append(f"Less Used Bonus[-{used_amount:.2f}]")

    # ---- NORMALIZED PLANET SCORES ----
    _nps_static_malefics = {'Sun', 'Mars', 'Saturn', 'Rahu', 'Ketu'}
    _nps_static_benefics = {'Jupiter', 'Venus', 'Mercury'}
    # (NATURAL_FRIENDSHIPS moved up)
    
    _nps_neecha_statuses = {'Neecham', 'Neechabhangam', 'Neechabhanga Raja Yoga'}
    _nps_moon_is_waxing = (paksha == 'Shukla') or (moon_phase_name == 'Purnima')

    nps_rows = []
    _nps_score_dict = {}  # planet -> raw final_ns value for use in Planet Strengths
    _suchama_score_dict = {}  # planet -> raw suchama value
    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        inv = phase5_data[p]['p5_inventory']
        p5_debt = phase5_data[p]['p5_current_debt']
        p_volume = phase5_data[p]['volume']
        p_capacity = capacity_dict.get(p, 100)

        total_good = sum(v for k, v in inv.items() if v > 0.001 and is_good_currency(k))
        total_bad  = sum(v for k, v in inv.items() if v > 0.001 and 'Bad' in k)
        # Jupiter Poison: treated as bad for NPS (poisonpenality = 3x)
        _nps_jp_poison = inv.get('Jupiter Poison', 0.0)
        if _nps_jp_poison > 0.001:
            total_good -= _nps_jp_poison
            total_bad += (poisonpenality - 1) * _nps_jp_poison
        self_bad   = inv.get(f'Bad {p}', 0.0)
        net_score  = total_good - total_bad

        _nps_st = planet_data[p].get('updated_status', '')
        if not _nps_st or _nps_st == '-':
            _nps_st = planet_data[p].get('status', '')
        is_neecha = _nps_st in _nps_neecha_statuses

        # Classify malefic/benefic
        if p == 'Moon':
            is_malefic = not _nps_moon_is_waxing
        elif p in _nps_static_malefics:
            is_malefic = True
        else:
            is_malefic = False

        # --- Determine case and calculate ---
        if p == 'Moon' and _nps_moon_is_waxing and not is_neecha:
            # Case A: Waxing Moon, NOT Negative Status
            abs_debt = abs(p5_debt)
            denom_val = total_good + abs_debt + total_bad
            if abs(denom_val) < 0.001:
                final_ns = 0.0
            else:
                final_ns = ((total_good - total_bad) / denom_val) * 100
            formula_type = f"CaseA: Waxing Moon [(Good {total_good:.2f} - Bad {total_bad:.2f}) / (Good {total_good:.2f} + |Debt| {abs_debt:.2f} + Bad {total_bad:.2f})] x100 = {final_ns:.2f}"

        elif p == 'Moon' and _nps_moon_is_waxing and is_neecha:
            # Case B: Waxing Moon, IS Negative Status
            swapped_debt = -1 * p5_debt
            denom_val = total_good + swapped_debt
            if abs(denom_val) < 0.001:
                final_ns = 0.0
            else:
                final_ns = ((total_good - swapped_debt) / denom_val) * 120
            formula_type = f"CaseB: [(TG{total_good:.2f}-SD{swapped_debt:.2f})/(TG{total_good:.2f}+SD{swapped_debt:.2f})]*120"

        elif not is_malefic and is_neecha:
            # Case C: Benefic, IS Negative Status
            denom_val = p_capacity * 1.2
            if abs(denom_val) < 0.001:
                final_ns = 0.0
            else:
                final_ns = (net_score / denom_val) * 120
            formula_type = f"CaseC: (Net{net_score:.2f}/(Cap{p_capacity}*1.2))*120"

        elif is_malefic and is_neecha:
            # Case D: Malefic, IS Negative Status
            denom_val = p_capacity * 1.2
            if abs(denom_val) < 0.001:
                final_ns = 0.0
            else:
                final_ns = ((net_score + self_bad) / denom_val) * 120
            formula_type = f"CaseD: ((Net{net_score:.2f}+SB{self_bad:.2f})/(Cap{p_capacity}*1.2))*120"

        elif not is_malefic and not is_neecha:
            # Case E: Benefic, NOT Negative Status
            if abs(p_volume) < 0.001:
                final_ns = 0.0
            else:
                final_ns = (net_score / p_volume) * 100
            formula_type = f"CaseE: (Net{net_score:.2f}/Vol{p_volume:.2f})*100"

        else:
            # Case F: Malefic, NOT Negative Status
            if abs(p_volume) < 0.001:
                final_ns = 0.0
            else:
                if p == 'Ketu':
                    final_ns = (net_score / p_volume) * 100
                    formula_type = f"CaseF: (Net{net_score:.2f}/Vol{p_volume:.2f})*100 [Ketu: SB excluded]"
                else:
                    final_ns = ((net_score + self_bad) / p_volume) * 100
                    formula_type = f"CaseF: ((Net{net_score:.2f}+SB{self_bad:.2f})/Vol{p_volume:.2f})*100"

        # KHS Calculation (Capped at 20) for NPS
        _khs_ruled = planet_ruled_signs.get(p, [])
        if not _khs_ruled:
            _khs_val = 0.0
        else:
            _khs_total = sum(aspect_score.get(rs, 0.0) + occupant_score.get(rs, 0.0) for rs in _khs_ruled)
            _khs_avg = _khs_total / len(_khs_ruled)
            # Old logic: _khs_val = min(_khs_avg / 10.0, 20.0)
            # New logic: Multiply by 2. No lower limit of 0, so negative values persist. 
            # Still apply cap of 20 on the positive side. 
            raw_khs = (_khs_avg / 10.0) * 2
            _khs_val = min(raw_khs, 20.0)

        _ns_without_khs = final_ns  # capture score before KHS

        final_ns += _khs_val
        if abs(_khs_val) > 0.001:
            formula_type += f" + KHS({_khs_val:.2f})"

        _nps_score_dict[p] = final_ns

        # Maraivu percentage lookup
        p_house = planet_house_map.get(p, 0)
        m_pct = maraivu_percentage.get(p, {}).get(p_house, None)
        m_pct_str = f"{m_pct}%" if m_pct is not None else "-"

        # Maraivu Adjusted Score and Suchama Score
        benefic_set = {'Moon', 'Mercury', 'Jupiter', 'Venus'}
        malefic_set = {'Sun', 'Mars', 'Saturn', 'Rahu', 'Ketu'}

        # Friendship groups for malefic maraivu logic
        group_a = {'Sun', 'Moon', 'Mars', 'Jupiter', 'Ketu'}
        group_b = {'Venus', 'Saturn', 'Mercury', 'Rahu'}

        suchama_str = "0"
        
        # New: Need to capture the Adjusted Score for use in House Points later
        final_adjusted_score = 0.0

        # Compute updated maraivu % (reduced by status) for Adjusted Score only
        _um_pct = m_pct  # default to raw maraivu %
        if m_pct is not None and m_pct > 0:
            _um_status = planet_status_map.get(p, '-')
            if _um_status == 'Uchcham':
                _um_pct = m_pct * 0.50
            elif _um_status == 'Moolathirigonam':
                _um_pct = m_pct * 0.60
            elif _um_status == 'Aatchi':
                _um_pct = m_pct * 0.70
            else:
                # Check Friend's House
                _um_lagna_sign = get_sign(lagna_sid)
                _um_house_sign = sign_names[(sign_names.index(_um_lagna_sign) + (p_house - 1)) % 12]
                _um_house_lord = sign_lords[sign_names.index(_um_house_sign)]
                _um_p_grp = 'A' if p in group_a else 'B'
                _um_l_grp = 'A' if _um_house_lord in group_a else 'B'
                if _um_p_grp == _um_l_grp:
                    _um_pct = m_pct * 0.75

        if p in benefic_set:
            # Benefic logic
            if _um_pct is not None:
                adjusted = (final_ns / 2.0) + (final_ns * (100 - _um_pct) / 100.0) / 2.0
            else:
                adjusted = final_ns
            final_adjusted_score = adjusted
            adjusted_str = f"{adjusted:.2f}"

        elif p in malefic_set:
            # Malefic logic
            p_sign = planet_sign_map.get(p, 'Aries')
            sthana_val = sthana_bala_dict.get(p, [0]*12)[sign_names.index(p_sign)]
            p_dig_bala = planet_data[p].get('dig_bala') or 0
            p_status = planet_status_map.get(p, '-')

            if m_pct is not None:
                # Check friendly house: house lord in same group as planet
                lagna_sign = get_sign(lagna_sid)
                house_sign = sign_names[(sign_names.index(lagna_sign) + (p_house - 1)) % 12]
                house_lord = sign_lords[sign_names.index(house_sign)]
                planet_group = 'A' if p in group_a else 'B'
                lord_group = 'A' if house_lord in group_a else 'B'
                is_friendly_house = (planet_group == lord_group)

                # Check if planet itself has positive status
                has_positive_status = p_status in ('Uchcham', 'Aatchi', 'Moolathirigonam')

                if is_friendly_house or has_positive_status:
                    # No reduction
                    adjusted = final_ns
                    # Suchama Score uses original m_pct (not reduced)
                    # Sun gets suchama only from digbala, not from maraivu
                    if p == 'Sun':
                        suchama = 0.0
                    else:
                        suchama = (sthana_val / 100.0) * m_pct
                else:
                    # Reduce using updated maraivu %
                    if final_ns < 0:
                         adjusted = (final_ns / 2.0) + (final_ns * (100 + _um_pct) / 100.0) / 2.0
                    else:
                         adjusted = (final_ns / 2.0) + (final_ns * (100 - _um_pct) / 100.0) / 2.0
                    suchama = 0.0

                # Step 1: If Digbala > 92%, add 0.5 * Sthana Balam
                if p_dig_bala > 92:
                    suchama += 0.5 * sthana_val

                # Step 2: If Saturn or Mars has Neecham, add 0.5 * Sthana Balam
                if p in ('Saturn', 'Mars') and p_status == 'Neecham':
                    suchama += 0.5 * sthana_val

                suchama_str = f"{suchama:.2f}"
                _suchama_score_dict[p] = suchama
                final_adjusted_score = adjusted
                adjusted_str = f"{adjusted:.2f}"
            else:
                # No maraivu detected
                adjusted = final_ns
                final_adjusted_score = adjusted
                adjusted_str = f"{adjusted:.2f}"
                suchama = 0.0

                # Step 1: If Digbala > 92%, add 0.5 * Sthana Balam
                if p_dig_bala > 92:
                    suchama += 0.5 * sthana_val

                # Step 2: If Saturn or Mars has Neecham, add 0.5 * Sthana Balam
                if p in ('Saturn', 'Mars') and p_status == 'Neecham':
                    suchama += 0.5 * sthana_val

                suchama_str = f"{suchama:.2f}"
                _suchama_score_dict[p] = suchama
        else:
            adjusted_str = "-"
            final_adjusted_score = final_ns # Default fallback if anything breaks 
            
        # Store adjusted score for House Points Calculation
        _nps_score_dict[p + '_adjusted'] = final_adjusted_score

        # ---- RAHU SCORE CALCULATION ----
        if p == 'Rahu':
            _rahu_total_bad_p5 = total_bad  # total Bad currency including self bad from phase 5
            _rahu_notes_parts = []

            # Step 1: Quantise Maraivu Adjusted Score to 80
            _rahu_base = (final_adjusted_score / 100.0) * 80
            _rahu_total = _rahu_base
            _rahu_notes_parts.append(f"Step1: ({final_adjusted_score:.2f}/100)*80={_rahu_base:.2f}")

            # Step 2: Adding score based on house lord status
            _rahu_sign = planet_sign_map.get('Rahu', 'Aries')
            _rahu_house_lord = get_sign_lord(_rahu_sign)
            _rahu_hl_status = planet_status_map.get(_rahu_house_lord, '-')
            if _rahu_hl_status == 'Uchcham':
                _rahu_total += 60
                _rahu_notes_parts.append(f"Step2: Lord {_rahu_house_lord} Uchcham +60")
            elif _rahu_hl_status == 'Moolathirigonam':
                _rahu_total += 48
                _rahu_notes_parts.append(f"Step2: Lord {_rahu_house_lord} Moolathirigonam +48")
            elif _rahu_hl_status == 'Aatchi':
                _rahu_total += 36
                _rahu_notes_parts.append(f"Step2: Lord {_rahu_house_lord} Aatchi +36")
            else:
                _rahu_notes_parts.append(f"Step2: Lord {_rahu_house_lord} status={_rahu_hl_status}, no bonus")

            # Step 3: Bonus if Rahu is in its favourite house
            _rahu_fav_signs = {'Aries', 'Taurus', 'Cancer', 'Virgo', 'Libra', 'Sagittarius', 'Capricorn', 'Pisces'}
            if _rahu_sign in _rahu_fav_signs:
                _rahu_total += 10
                _rahu_notes_parts.append(f"Step3: Fav house {_rahu_sign} +10")
            else:
                _rahu_notes_parts.append(f"Step3: {_rahu_sign} not fav, no bonus")

            # Step 4: Bonus if Rahu is in friend's house for ascendant lagna
            _rahu_friend_list = {
                'Aries': {'Cancer', 'Leo', 'Scorpio', 'Sagittarius', 'Pisces'},
                'Taurus': {'Gemini', 'Virgo', 'Libra', 'Capricorn', 'Aquarius'},
                'Gemini': {'Taurus', 'Virgo', 'Libra', 'Capricorn', 'Aquarius'},
                'Cancer': {'Aries', 'Leo', 'Scorpio', 'Sagittarius', 'Pisces'},
                'Leo': {'Aries', 'Cancer', 'Scorpio', 'Sagittarius', 'Pisces'},
                'Virgo': {'Taurus', 'Gemini', 'Libra', 'Capricorn', 'Aquarius'},
                'Libra': {'Taurus', 'Gemini', 'Virgo', 'Capricorn', 'Aquarius'},
                'Scorpio': {'Aries', 'Cancer', 'Leo', 'Sagittarius', 'Pisces'},
                'Sagittarius': {'Aries', 'Cancer', 'Leo', 'Scorpio', 'Pisces'},
                'Capricorn': {'Taurus', 'Gemini', 'Virgo', 'Libra', 'Aquarius'},
                'Aquarius': {'Taurus', 'Gemini', 'Virgo', 'Libra', 'Capricorn'},
                'Pisces': {'Aries', 'Cancer', 'Leo', 'Scorpio', 'Sagittarius'},
            }
            _rahu_asc_sign = get_sign(lagna_sid)
            _rahu_friends = _rahu_friend_list.get(_rahu_asc_sign, set())
            if _rahu_sign in _rahu_friends:
                _rahu_total += 30
                _rahu_notes_parts.append(f"Step4: {_rahu_sign} friend of lagna {_rahu_asc_sign} +30")
            else:
                _rahu_notes_parts.append(f"Step4: {_rahu_sign} not friend of lagna {_rahu_asc_sign}, no bonus")

            rahu_score_str = f"{_rahu_total:.2f}"
            rahu_notes_str = " | ".join(_rahu_notes_parts)
        else:
            rahu_score_str = "-"
            rahu_notes_str = "-"

        # ---- PLANET HAPPINESS SCORE (Naisargika Graha Maitri) ----
        _hap_score = 0
        _hap_notes = []
        # A. Residence Score
        _hap_planet_sign = planet_sign_map.get(p, 'Aries')
        _hap_sign_lord = get_sign_lord(_hap_planet_sign)
        _hap_own_signs = planet_ruled_signs.get(p, [])
        if _hap_planet_sign in _hap_own_signs:
            _hap_score += 2
            _hap_notes.append(f"In Own({_hap_planet_sign}) [+2]")
        else:
            _hap_rel = check_friendship(p, _hap_sign_lord)
            if _hap_rel == 'Friend':
                _hap_score += 2
                _hap_notes.append(f"In {_hap_sign_lord}(Fr) House [+2]")
            elif _hap_rel == 'Enemy':
                _hap_score -= 2
                _hap_notes.append(f"In {_hap_sign_lord}(En) House [-2]")
            else:
                _hap_notes.append(f"In {_hap_sign_lord}(Ne) House [+0]")
        # B. Aspect Score (Incoming Clones within 22 degrees)
        _hap_planet_L = phase5_data[p]['L']
        for _hap_cl in all_initial_clones:
            _hap_source = _hap_cl['parent']
            if _hap_source == p:
                continue  # skip own clones
            _hap_cl_L = _hap_cl['L']
            _hap_diff = abs(_hap_planet_L - _hap_cl_L)
            if _hap_diff > 180:
                _hap_diff = 360 - _hap_diff
            if _hap_diff <= 22:
                _hap_asp_rel = check_friendship(p, _hap_source)
                if _hap_asp_rel == 'Friend':
                    _hap_score += 1
                    _hap_notes.append(f"Asp by {_hap_source}(Fr) [+1]")
                elif _hap_asp_rel == 'Enemy':
                    _hap_score -= 1
                    _hap_notes.append(f"Asp by {_hap_source}(En) [-1]")
                else:
                    _hap_notes.append(f"Asp by {_hap_source}(Ne) [+0]")
        _hap_score_str = str(_hap_score)
        _hap_notes_str = " | ".join(_hap_notes) if _hap_notes else "-"

        # Currency % breakdown from Phase 5 inventory
        # Waxing Moon: Currency / (Good + Bad + |Debt|) * 100
        # Others: Currency / volume * 100
        _cur_pct_parts = []
        if p == 'Moon' and _nps_moon_is_waxing:
            _cur_denom = total_good + total_bad + abs(p5_debt)
            _cur_denom = _cur_denom if _cur_denom > 0.001 else 1.0
        else:
            _cur_denom = p_volume if p_volume > 0.001 else 1.0
        for _ck in sorted(inv.keys(), key=lambda x: abs(inv[x]), reverse=True):
            _cv = inv[_ck]
            if abs(_cv) > 0.001:
                _cpct = (abs(_cv) / _cur_denom) * 100.0
                _cur_pct_parts.append(f"{_ck}({_cpct:.1f}%)")
        _cur_pct_str = ", ".join(_cur_pct_parts) if _cur_pct_parts else "-"

        nps_rows.append([p, f"{net_score:.2f}", f"{self_bad:.2f}", formula_type, f"{final_ns:.2f}", f"{_ns_without_khs:.2f}", _cur_pct_str, m_pct_str, adjusted_str, suchama_str, rahu_score_str, rahu_notes_str, _hap_score_str, _hap_notes_str])

    df_normalized_planet_scores = pd.DataFrame(nps_rows,
        columns=['Planet', 'Net Score', 'Self Bad', 'Formula Type', 'Final Normalized Score', 'Normalised without KHS', 'Currency %', 'Maraivu %', 'Maraivu Adjusted Score', 'Suchama Score', 'Rahu Score', 'Rahu Notes', 'Happiness Score', 'Happiness Notes'])

    # ---- NEECHAM STATUS UPGRADE BASED ON FINAL NORMALIZED SCORE ----
    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        _p_status = planet_data[p].get('status', '-')
        _p_updated = planet_data[p].get('updated_status', '-')
        if _p_status == 'Neecham':
            _p_fns = _nps_score_dict.get(p, 0.0)
            if _p_updated == 'Neechabhangam':
                # Already Neechabhangam — promote to Raja Yoga if score > 115
                if _p_fns > 115:
                    planet_data[p]['updated_status'] = 'Neechabhanga Raja Yoga'
                    for row in rows:
                        if row[0] == p:
                            row[11] = 'Neechabhanga Raja Yoga'
            elif _p_updated in ('-', '', None):
                # No updated status yet — assign based on score
                if _p_fns > 115:
                    planet_data[p]['updated_status'] = 'Neechabhanga Raja Yoga'
                    for row in rows:
                        if row[0] == p:
                            row[11] = 'Neechabhanga Raja Yoga'
                elif _p_fns > 75:
                    planet_data[p]['updated_status'] = 'Neechabhangam'
                    for row in rows:
                        if row[0] == p:
                            row[11] = 'Neechabhangam'

    # ---- Apply Rahu Score directly as occupant score for Rahu's house ----
    _rahu_occ_sign = planet_sign_map.get('Rahu', 'Aries')
    _rahu_bad_penalty = _rahu_total_bad_p5 * 0.80
    _rahu_occupant_val = _rahu_total - _rahu_bad_penalty
    occupant_score[_rahu_occ_sign] += _rahu_occupant_val
    occupant_notes[_rahu_occ_sign].append(f"Rahu(RahuScore={_rahu_total:.2f} - BadPenalty={_rahu_bad_penalty:.2f} => {_rahu_occupant_val:.2f})")

    # ---- PLANET STRENGTHS ANALYSIS ----
    planet_strength_rows = []
    planet_final_strengths = {}
    _ps_planets = ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn']

    def _ps_own_asp(p_name):
        ruled = planet_ruled_signs.get(p_name, [])
        if not ruled:
            return 0.0
        p_L = phase5_data[p_name]['L']
        for cl in all_leftover_clones:
            if cl['parent'] == p_name:
                t_sign = get_sign((p_L + (cl['offset'] - 1) * 30) % 360)
                if t_sign in ruled:
                    return 10.0
        return 0.0

    _overridden_sthana = {}   # track planets whose sthana was overridden
    _planet_maraivu_adj_strengths = {} 

    for _ps_p in _ps_planets:
        _db = planet_data[_ps_p].get('dig_bala') or 0
        _sb = planet_data[_ps_p].get('sthana') or 0

        # Override Sthana with normalized score for negative-status planets
        _neg_statuses = {'Neecham', 'Neechabhangam', 'Neechabhanga Raja Yoga'}
        _ps_st = planet_data[_ps_p].get('status', '')
        _ps_ust = planet_data[_ps_p].get('updated_status', '')
        _is_negative = _ps_st in _neg_statuses or _ps_ust in _neg_statuses
        if _is_negative:
            # Use Final Normalized Score from normalized planet scores (scale 0-120)
            _nps_val = _nps_score_dict.get(_ps_p, 0.0)
            _sb = _nps_val  # will be converted to 40 or 60 scale below via (val/120)*weight
            _overridden_sthana[_ps_p] = _sb

        # Parivardhana Yoga: override Sthana with swapped sign score
        _pari_note = ''
        _pari_entry = parivardhana_map.get(_ps_p, '')
        if _pari_entry and not _is_negative:
            _partner_name = _pari_entry.split(' (')[0]
            _partner_sign = planet_sign_map.get(_partner_name, '')
            if _partner_sign and _partner_sign in sign_names:
                _partner_idx = sign_names.index(_partner_sign)
                _sb = sthana_bala_dict.get(_ps_p, [0]*12)[_partner_idx]
                _pari_note = f'(Pari->{_partner_sign})'

        _asp_val = _ps_own_asp(_ps_p)
        _rh = phase5_data[_ps_p]['rasi_house']

        # House Lord Status bonus/penalty
        _ps_sign = phase5_data[_ps_p]['sign']
        _ps_lord = get_sign_lord(_ps_sign)
        _lord_st = planet_status_map.get(_ps_lord, '-') if _ps_lord else '-'
        if _lord_st == 'Uchcham':
            _hl_adj = 20.0
        elif _lord_st == 'Moolathirigonam':
            _hl_adj = 16.0
        elif _lord_st == 'Aatchi':
            _hl_adj = 12.0
        elif _lord_st == 'Neecham':
            _hl_adj = -20.0
        else:
            _hl_adj = 0.0

        # If the planet itself is Neecha, do not award HLord bonus (ucham/moola/aatchi)
        if _is_negative and _hl_adj > 0:
            _hl_adj = 0.0

        if _hp_is_malefic(_ps_p):
            # Malefic: Dig 60%, Sthana 40%, Asp 10%, Kendra 5%
            s_dig = (_db / 100.0) * 60.0
            if _is_negative:
                s_sth = (_sb / 100.0) * 40.0
            else:
                s_sth = (_sb / 100.0) * 40.0
            base_total = s_dig + s_sth + _asp_val + _hl_adj
            s_bonus = 5.0 if _rh in (1, 4, 7, 10) else 0.0
        else:
            # Benefic: Dig 40%, Sthana 60%, Asp 10%, Kona 5%
            s_dig = (_db / 100.0) * 40.0
            if _is_negative:
                s_sth = (_sb / 100.0) * 60.0
            else:
                s_sth = (_sb / 100.0) * 60.0
            base_total = s_dig + s_sth + _asp_val + _hl_adj
            s_bonus = 5.0 if _rh in (1, 5, 9) else 0.0

        final = base_total + s_bonus
        planet_final_strengths[_ps_p] = final
        brkdn = (f"Dig:{s_dig:.2f} + Sthana:{s_sth:.2f}{_pari_note} + "
                 f"Asp:{_asp_val:.2f} + "
                 f"HLord:{_hl_adj:+.2f}({_ps_lord}={_lord_st}) + Bonus:{s_bonus:.2f}")

        # ── Updated Maraivu Percentage & Adjusted Strength ──
        _ps_rh = phase5_data[_ps_p]['rasi_house']
        _ps_base_maraivu = maraivu_percentage.get(_ps_p, {}).get(_ps_rh, 0)

        _ps_updated_maraivu = _ps_base_maraivu  # start with base
        if _ps_base_maraivu > 0:
            # Check status-based reduction (reduce BY x% of the base value)
            _ps_planet_status = planet_status_map.get(_ps_p, '-')
            if _ps_planet_status == 'Uchcham':
                _ps_updated_maraivu = _ps_base_maraivu * 0.50       # reduce by 50%
            elif _ps_planet_status == 'Moolathirigonam':
                _ps_updated_maraivu = _ps_base_maraivu * 0.60       # reduce by 40%
            elif _ps_planet_status == 'Aatchi':
                _ps_updated_maraivu = _ps_base_maraivu * 0.70       # reduce by 30%
            else:
                # Check Friend's House: planet and house lord in same group
                _ps_fr_group_a = {'Sun', 'Moon', 'Mars', 'Jupiter', 'Ketu'}
                _ps_fr_group_b = {'Venus', 'Mercury', 'Saturn', 'Rahu'}
                _ps_house_sign = sign_names[(sign_names.index(get_sign(lagna_sid)) + (_ps_rh - 1)) % 12]
                _ps_house_lord = sign_lords[sign_names.index(_ps_house_sign)]
                _ps_p_grp = 'A' if _ps_p in _ps_fr_group_a else 'B'
                _ps_l_grp = 'A' if _ps_house_lord in _ps_fr_group_a else 'B'
                if _ps_p_grp == _ps_l_grp:
                    _ps_updated_maraivu = _ps_base_maraivu * 0.75   # reduce by 25%

        # Maraivu Adjusted Strength
        _ps_adj_strength = (final / 2.0) + (final * (100 - _ps_updated_maraivu) / 200.0)
        _planet_maraivu_adj_strengths[_ps_p] = _ps_adj_strength

        planet_strength_rows.append([_ps_p, f"{final:.2f}", f"{_ps_adj_strength:.2f}", brkdn])

    df_planet_strengths = pd.DataFrame(planet_strength_rows,
        columns=['Planet', 'Total Strength', 'Maraivu Adjusted Strength', 'Score Breakdown'])

    # Update planet_data and rows with overridden Sthana Bala for negative-status planets
    # Update planet_data and rows with overridden Sthana Bala for negative-status planets
    if _overridden_sthana:
        for row in rows:
            p_name = row[0]
            if p_name in _overridden_sthana:
                row[9] = f"{_overridden_sthana[p_name]:.2f}%"
                planet_data[p_name]['sthana'] = _overridden_sthana[p_name]
        # Recreate df_planets so it reflects the overridden Sthana Bala values
        df_planets = pd.DataFrame(rows, columns=['Planet','Deg','Sign','Nakshatra','Pada','Ld/SL','Vargothuva',
                                                 'Parivardhana',
                                                 'Dig Bala (%)','Sthana Bala (%)','Status','Updated Status',
                                                 'Volume', 'Default Currencies', 'Debt'])
    # ---- END PLANET STRENGTHS ----

    # ── 4. BUILD DATAFRAME ──
    hp_rows = []
    _house_total_points = {}  # house_number -> raw total_hp
    _house_planetary_scores = {}  # house_number -> house_planetary_score
    _hp_lagna_idx = sign_names.index(get_sign(lagna_sid))
    for h_num in range(1, 13):
        s = sign_names[(_hp_lagna_idx + h_num - 1) % 12]
        a_src = ', '.join(aspect_sources[s]) if aspect_sources[s] else '-'
        o_src = ', '.join(occupant_notes[s]) if occupant_notes[s] else '-'

        house_planetary_score = aspect_score[s] + occupant_score[s]

        # House Lord Score = (lord strength / 2) + (Maraivu Adjusted Score / 2)
        lord = get_sign_lord(s)
        lord_strength = _planet_maraivu_adj_strengths.get(lord, 0.0)
        # Old: lord_norm_score = _nps_score_dict.get(lord, 0.0)
        # New: Use Maraivu adjusted score
        lord_norm_score = _nps_score_dict.get(lord + '_adjusted', 0.0)
        
        hl_score = (lord_strength / 2.0) + (lord_norm_score / 2.0)
        hl_notes = f"{lord}: Str({lord_strength/2.0:.2f}) + AdjNormScore({lord_norm_score/2.0:.2f})"

        # Total House Points = (House Planetary Score / 2) + (House Lord Score / 2)
        total_hp = (house_planetary_score / 2.0) + (hl_score / 2.0)
        total_hp_notes = f"HPS({house_planetary_score/2.0:.2f}) + HLS({hl_score/2.0:.2f})"

        _house_total_points[h_num] = total_hp
        _house_planetary_scores[h_num] = house_planetary_score

        # ---- HOUSE HAPPINESS SCORE (Naisargika Graha Maitri) ----
        _hh_score = 0
        _hh_notes = []
        _hh_lord = get_sign_lord(s)
        # B. Occupant Score
        _hh_occupants = [occ for occ in house_planets_rasi[h_num] if occ not in ('Asc',)]
        for _hh_occ in _hh_occupants:
            if _hh_occ == _hh_lord:
                continue  # lord in own house, skip as occupant
            _hh_occ_rel = check_friendship(_hh_lord, _hh_occ)
            if _hh_occ_rel == 'Friend':
                _hh_score += 1
                _hh_notes.append(f"Occ {_hh_occ}(Fr) [+1]")
            elif _hh_occ_rel == 'Enemy':
                _hh_score -= 1
                _hh_notes.append(f"Occ {_hh_occ}(En) [-1]")
            else:
                _hh_notes.append(f"Occ {_hh_occ}(Ne) [+0]")
        # C. Aspect Score (Clones falling into this house sign)
        for _hh_cl in all_initial_clones:
            _hh_cl_sign = get_sign(_hh_cl['L'])
            if _hh_cl_sign == s:
                _hh_asp_source = _hh_cl['parent']
                _hh_asp_rel = check_friendship(_hh_lord, _hh_asp_source)
                if _hh_asp_rel == 'Friend':
                    _hh_score += 1
                    _hh_notes.append(f"Asp {_hh_asp_source}(Fr) [+1]")
                elif _hh_asp_rel == 'Enemy':
                    _hh_score -= 1
                    _hh_notes.append(f"Asp {_hh_asp_source}(En) [-1]")
                else:
                    _hh_notes.append(f"Asp {_hh_asp_source}(Ne) [+0]")
        _hh_score_str = str(_hh_score)
        _hh_notes_str = " | ".join(_hh_notes) if _hh_notes else "-"

        hp_rows.append([h_num, s, f"{aspect_score[s]:.2f}", a_src, f"{occupant_score[s]:.2f}", o_src,
                        f"{house_planetary_score:.2f}",
                        f"{hl_score:.2f}", hl_notes,
                        f"{total_hp:.2f}", total_hp_notes, _hh_score_str, _hh_notes_str])

    df_house_points = pd.DataFrame(hp_rows,
        columns=['House', 'House Sign', 'Aspect Score', 'Aspect Sources', 'Occupant Score', 'Occupant Notes',
                 'House Planetary Score',
                 'House Lord Score', 'House Lord Score Notes',
                 'Total House Points', 'Total House Points Notes', 'Happiness Score', 'Happiness Notes'])
    # ---- END HOUSE POINTS ----

    df_planets = pd.DataFrame(rows, columns=['Planet','Deg','Sign','Nakshatra','Pada','Ld/SL','Vargothuva',
                                             'Parivardhana',
                                             'Dig Bala (%)','Sthana Bala (%)','Status','Updated Status',
                                             'Volume', 'Default Currencies', 'Debt'])

    df_rasi = pd.DataFrame([[f"House {h}", get_sign((lagna_sid+(h-1)*30)%360), 
                             ', '.join(sorted(house_planets_rasi[h])) if house_planets_rasi[h] else 'Empty'] 
                            for h in range(1,13)], columns=['House','Sign','Planets'])

    house_planets_nav = defaultdict(list)
    # Add Ascendant to Navamsa House 1
    house_planets_nav[1].append("Asc")
    for p,L in lon_sid.items():
        nav_lon = (L*9) % 360
        nav_h = (int(nav_lon/30) - int(nav_lagna/30)) % 12 + 1
        house_planets_nav[nav_h].append(p.capitalize())
    
    df_nav = pd.DataFrame([[f"House {h}", get_sign((nav_lagna+(h-1)*30)%360), 
                            ', '.join(sorted(house_planets_nav[h])) if house_planets_nav[h] else 'Empty'] 
                           for h in range(1,13)], columns=['House','Sign','Planets'])

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

    moon_lon = lon_sid['moon']
    idx, bal = generate_vimshottari_dasa(moon_lon)
    full_first = years[idx]; passed = full_first - bal
    dasa_start = utc_dt - timedelta(days=passed*365.25)
    dasa = generate_periods(dasa_start, idx, 120, 'dasa', max_depth)
    dasa_filtered = filter_from_birth(dasa, utc_dt)

    depth_map = {1:'Dasa only',2:'Dasa + Bhukti',3:'Dasa + Bhukti + Anthara',
                 4:'Dasa + Bhukti + Anthara + Sukshma',5:'Dasa + Bhukti + Anthara + Sukshma + Prana',
                 6:'Dasa + Bhukti + Anthara + Sukshma + Prana + Sub-Prana'}

    # ====== LAGNA POINT SCORE SIMULATION ======
    # 1. Preparation: Build "Universe of Pots"
    _sim_malefic_names = {'Saturn', 'Mars', 'Sun', 'Rahu', 'Ketu'}

    def _sim_is_malefic(pname):
        if pname in _sim_malefic_names:
            return True
        if pname == 'Moon':
            return phase5_data['Moon']['bad_inv'] > 0.001
        return False

    universe_pots = []

    # A. Real Planet Pots
    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        pot_inv = copy.deepcopy(dict(phase5_data[p]['p5_inventory']))
        pot_vol = phase5_data[p]['volume']
        pot_L = phase5_data[p]['L']
        universe_pots.append({
            'name': p,
            'L': pot_L,
            'inventory': pot_inv,
            'volume': pot_vol,
            'is_malefic': _sim_is_malefic(p),
            'kind': 'real'
        })

    # B. Virtual Clone Pots (Leftover Aspects)
    for cl in all_leftover_clones:
        cl_inv = copy.deepcopy(dict(cl['inventory']))
        cl_vol = sum(v for v in cl_inv.values() if v > 0.001)
        parent = cl['parent']
        universe_pots.append({
            'name': f"Clone({parent}_{cl['offset']})",
            'L': cl['L'],
            'inventory': cl_inv,
            'volume': cl_vol,
            'is_malefic': _sim_is_malefic(parent),
            'kind': 'clone'
        })

    # 2. The Simulation: Lagna Pulling
    sim_lagna_L = lagna_sid
    sim_debt = -100.0
    sim_gained_inv = defaultdict(float)
    sim_sources = {}  # currency_key -> list of "PotName(amount)" strings
    sim_good_from_malefic = defaultdict(float)  # good currency amounts sourced from malefic pots

    malefic_pots = [p for p in universe_pots if p['is_malefic']]
    benefic_pots = [p for p in universe_pots if not p['is_malefic']]
    ordered_pots = malefic_pots + benefic_pots

    for pot in ordered_pots:
        if sim_debt >= -0.001:
            break

        # Calculate gap
        raw_diff = abs(sim_lagna_L - pot['L'])
        if raw_diff > 180:
            raw_diff = 360 - raw_diff
        gap = int(raw_diff)

        if gap > 22:
            continue

        cap_pct = mix_dict.get(gap, 0)
        max_pull = pot['volume'] * (cap_pct / 100.0)
        remaining_cap = max_pull

        # Sort inventory: prioritise Good currency (higher rank first)
        sorted_currencies = sorted(
            [(k, v) for k, v in pot['inventory'].items() if v > 0.001],
            key=lambda x: get_p5_currency_rank_score(x[0]),
            reverse=True
        )

        for c_key, c_avail in sorted_currencies:
            if sim_debt >= -0.001 or remaining_cap <= 0.001:
                break
            needed = abs(sim_debt)
            take = min(needed, c_avail, remaining_cap)
            if take > 0.001:
                sim_gained_inv[c_key] += take
                sim_debt += take
                remaining_cap -= take
                # Track source
                src_label = pot['name']
                sim_sources.setdefault(c_key, []).append(f"{src_label}({take:.2f})")
                # Track good currency from malefic pots
                if pot['is_malefic'] and is_good_currency(c_key) and c_key != 'Jupiter Poison':
                    sim_good_from_malefic[c_key] += take

    # Post-sim: halve good currency that came from malefic pots
    for c_key, malefic_amount in sim_good_from_malefic.items():
        penalty = malefic_amount * 0.50
        sim_gained_inv[c_key] -= penalty

    # 3. Output Generation
    sim_good_total = sum(v for k, v in sim_gained_inv.items() if is_good_currency(k))
    sim_bad_total = sum(v for k, v in sim_gained_inv.items() if 'Bad' in k)
    # Jupiter Poison: treat as bad in Lagna simulation
    _sim_jp_poison = sim_gained_inv.get('Jupiter Poison', 0.0)
    if _sim_jp_poison > 0.001:
        sim_good_total -= _sim_jp_poison
        sim_bad_total += _sim_jp_poison
    sim_net_score = sim_good_total - sim_bad_total

    # Currency breakdown string
    breakdown_parts = []
    for k in sorted(sim_gained_inv.keys(), key=lambda x: get_p5_currency_rank_score(x), reverse=True):
        v = sim_gained_inv[k]
        if v > 0.001:
            breakdown_parts.append(f"{k}[{v:.2f}]")
    breakdown_str = ", ".join(breakdown_parts) if breakdown_parts else "-"

    # Notes: source details per currency
    notes_parts = []
    for k in sorted(sim_sources.keys(), key=lambda x: get_p5_currency_rank_score(x), reverse=True):
        entries = sim_sources[k]
        notes_parts.append(f"{k} from " + ", ".join(entries))
    notes_str = "; ".join(notes_parts) if notes_parts else "-"

    remaining_debt_str = f"{sim_debt:.2f}" if abs(sim_debt) >= 0.01 else '0.00'

    # ====== END LAGNA POINT SCORE SIMULATION ======

    # ====== NAVAMSA LAGNA SCORE SIMULATION ======
    # Simulates an imaginary planet at the Navamsa Lagna with -100 debt,
    # pulling currency from Navamsa Phase 3 ecosystem.

    # Step A: Setup & Cloning
    _nav_lagna_house = 1  # Navamsa Lagna is always House 1
    nav_sim_data = {}
    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        nav_sim_data[p] = {
            'inventory': copy.deepcopy(dict(navamsa_phase3_data[p]['navp3_inventory'])),
            'debt': navamsa_phase3_data[p]['navp3_current_debt'],
            'volume': navamsa_phase3_data[p]['nav_volume'],
            'house': navamsa_phase3_data[p]['nav_house'],
        }
    nav_sim_house_pot = copy.deepcopy(dict(house_pot))

    sim_nav_debt = -100.0
    sim_nav_gained = defaultdict(float)
    sim_nav_sources = {}

    # Determine malefic/benefic for navamsa sim
    _nav_sim_malefic_names = {'Saturn', 'Mars', 'Sun', 'Rahu', 'Ketu'}
    _nav_moon_bad = nav_sim_data['Moon']['inventory'].get('Bad Moon', 0.0)
    _nav_moon_is_malefic = _nav_moon_bad > 0.001
    _nav_moon_is_waxing = (paksha == 'Shukla') or (moon_phase_name == 'Purnima')

    def _nav_sim_is_malefic(pname):
        if pname in _nav_sim_malefic_names:
            return True
        if pname == 'Moon':
            return _nav_moon_is_malefic
        return False

    def _nav_sim_is_benefic(pname):
        if pname in ('Jupiter', 'Venus', 'Mercury'):
            return True
        if pname == 'Moon' and not _nav_moon_is_malefic:
            return True
        return False

    # Step B: Interaction with House Pot for Navamsa Lagna House
    _nav_pot_val = nav_sim_house_pot.get(_nav_lagna_house, 0.0)
    if _nav_pot_val > 0.001 and sim_nav_debt < -0.001:
        _nav_pot_take = min(abs(sim_nav_debt), _nav_pot_val)
        sim_nav_debt += _nav_pot_take
        sim_nav_gained['Good Moon'] += _nav_pot_take
        nav_sim_house_pot[_nav_lagna_house] -= _nav_pot_take
        sim_nav_sources.setdefault('Good Moon', []).append(f"HousePot_H{_nav_lagna_house}({_nav_pot_take:.2f})")

    # Step C: Interaction with Planets (Sucking Phase)
    # Iterate through all planets; only those in the same Navamsa house as Lagna (house 1)
    _nav_planet_order = ['Saturn', 'Rahu', 'Sun', 'Mars', 'Ketu', 'Moon', 'Jupiter', 'Venus', 'Mercury']

    for _nav_tp in _nav_planet_order:
        if sim_nav_debt >= -0.001:
            break

        _nav_tp_data = nav_sim_data[_nav_tp]

        # Location Check: must be in the same house as Navamsa Lagna
        if _nav_tp_data['house'] != _nav_lagna_house:
            continue

        # Benefic Debt Check: if benefic and its debt is deeper than lagna's remaining debt, skip
        if _nav_sim_is_benefic(_nav_tp):
            _nav_tp_abs_debt = abs(_nav_tp_data['debt']) if _nav_tp_data['debt'] < -0.001 else 0.0
            _nav_sim_abs_debt = abs(sim_nav_debt)
            if _nav_tp_abs_debt >= _nav_sim_abs_debt:
                continue

        # Currency Rules: malefic -> only bad currencies; benefic -> only good currencies
        _nav_is_mal = _nav_sim_is_malefic(_nav_tp)
        _nav_inv = _nav_tp_data['inventory']

        _nav_allowable = []
        for _nk, _nv in _nav_inv.items():
            if _nv <= 0.001:
                continue
            if _nav_is_mal:
                if 'Bad' in _nk:
                    _nav_allowable.append((_nk, _nv))
            else:
                if is_good_currency(_nk):
                    _nav_allowable.append((_nk, _nv))

        # Sort by rank (higher first)
        _nav_allowable.sort(key=lambda x: get_p5_currency_rank_score(x[0]), reverse=True)

        # Pull currency
        for _nc_key, _nc_avail in _nav_allowable:
            if sim_nav_debt >= -0.001:
                break
            _nav_needed = abs(sim_nav_debt)
            _nav_take = min(_nav_needed, _nc_avail)
            if _nav_take > 0.001:
                sim_nav_gained[_nc_key] += _nav_take
                sim_nav_debt += _nav_take
                _nav_inv[_nc_key] -= _nav_take
                sim_nav_sources.setdefault(_nc_key, []).append(f"{_nav_tp}({_nav_take:.2f})")

    # Output Calculation
    sim_nav_good_total = sum(v for k, v in sim_nav_gained.items() if is_good_currency(k))
    sim_nav_bad_total = sum(v for k, v in sim_nav_gained.items() if 'Bad' in k)
    sim_nav_net_score = sim_nav_good_total - sim_nav_bad_total

    # Currency breakdown string
    _nav_bd_parts = []
    for k in sorted(sim_nav_gained.keys(), key=lambda x: get_p5_currency_rank_score(x), reverse=True):
        v = sim_nav_gained[k]
        if v > 0.001:
            _nav_bd_parts.append(f"{k}[{v:.2f}]")
    nav_breakdown_str = ", ".join(_nav_bd_parts) if _nav_bd_parts else "-"

    # Notes: source details
    _nav_notes_parts = []
    for k in sorted(sim_nav_sources.keys(), key=lambda x: get_p5_currency_rank_score(x), reverse=True):
        entries = sim_nav_sources[k]
        _nav_notes_parts.append(f"{k} from " + ", ".join(entries))
    nav_notes_str = "; ".join(_nav_notes_parts) if _nav_notes_parts else "-"

    nav_remaining_debt_str = f"{sim_nav_debt:.2f}" if abs(sim_nav_debt) >= 0.01 else '0.00'

    # ====== END NAVAMSA LAGNA SCORE SIMULATION ======

    df_bonus = pd.DataFrame(
        [
            ['Lagna Score', f"{sim_net_score:.2f}", '-100.00', remaining_debt_str, breakdown_str, notes_str],
            ['Navamsa Lagna Score', f"{sim_nav_net_score:.2f}", '-100.00', nav_remaining_debt_str, nav_breakdown_str, nav_notes_str]
        ],
        columns=['Simulation', 'Score', 'Initial Debt', 'Remaining Debt', 'Currency Breakdown', 'Notes']
    )
    # ====== END LAGNA POINT SCORE SIMULATION ======

    # ====== LAGNA ANALYSIS TABLE ======
    _la_lagna_sign = get_sign(lagna_sid)
    _la_lagna_lord = get_sign_lord(_la_lagna_sign)

    # 1. Moon's Light
    _la_moon_inv = phase5_data['Moon']['p5_inventory']
    _la_moon_good = sum(v for k, v in _la_moon_inv.items() if v > 0.001 and is_good_currency(k))
    _la_moon_bad = sum(v for k, v in _la_moon_inv.items() if v > 0.001 and 'Bad' in k)
    # Jupiter Poison: treat as bad for Moon's Light
    _la_moon_jp_poison = _la_moon_inv.get('Jupiter Poison', 0.0)
    if _la_moon_jp_poison > 0.001:
        _la_moon_good -= _la_moon_jp_poison
        _la_moon_bad += _la_moon_jp_poison
    _la_moon_debt = phase5_data['Moon']['p5_current_debt']
    _la_moon_abs_debt = abs(_la_moon_debt) if _la_moon_debt < -0.001 else 0.0
    _la_moon_is_waxing = (paksha == 'Shukla') or (moon_phase_name == 'Purnima')
    if _la_moon_is_waxing:
        # Waxing: [(Total Good - Total Bad) / (Total Good + |Debt|)] x 100
        _la_moon_denom = _la_moon_good + _la_moon_abs_debt
        if abs(_la_moon_denom) > 0.001:
            _la_moon_score = ((_la_moon_good - _la_moon_bad) / _la_moon_denom) * 100.0
        else:
            _la_moon_score = 0.0
        _la_moon_notes = f"Waxing Moon [(Good {_la_moon_good:.2f} - Bad {_la_moon_bad:.2f}) / (Good {_la_moon_good:.2f} + |Debt| {_la_moon_abs_debt:.2f})] x100 = {_la_moon_score:.2f}"
    else:
        # Waning: [(Total Good - |Debt|) / (Total Good + |Debt|)] x 100
        _la_moon_total = _la_moon_good + _la_moon_abs_debt
        if _la_moon_total > 0.001:
            _la_moon_score = ((_la_moon_good - _la_moon_abs_debt) / _la_moon_total) * 100.0
        else:
            _la_moon_score = 0.0
        _la_moon_notes = f"Waning Moon [(Good {_la_moon_good:.2f} - Debt {_la_moon_abs_debt:.2f}) / (Good {_la_moon_good:.2f} + Debt {_la_moon_abs_debt:.2f})] x100 = {_la_moon_score:.2f}"

    # 2. Lagna Lord Maraivu Adj Score from NPS
    _la_ll_adj = _nps_score_dict.get(_la_lagna_lord + '_adjusted', 0.0)
    _la_ll_score = _la_ll_adj
    _la_ll_notes = f"{_la_lagna_lord} maraivu adj NPS = {_la_ll_adj:.2f}"

    # 3. Lagna Lord Strength (maraivu adj from Planet Strengths)
    _la_ll_str_raw = _planet_maraivu_adj_strengths.get(_la_lagna_lord, 0.0)
    _la_ll_str_score = _la_ll_str_raw
    _la_ll_str_notes = f"{_la_lagna_lord} maraivu adj strength = {_la_ll_str_raw:.2f}"

    # 4. Lagna Lord Shukshama Strength
    _la_ll_suchama = _suchama_score_dict.get(_la_lagna_lord, 0.0)
    _la_ll_suchama_score = _la_ll_suchama
    _la_ll_suchama_notes = f"{_la_lagna_lord} suchama = {_la_ll_suchama:.2f}"

    # 5. 1st House Planetary Score
    _la_h1_raw = _house_planetary_scores.get(1, 0.0)
    _la_h1_score = _la_h1_raw
    _la_h1_notes = f"House 1 Planetary Score = {_la_h1_raw:.2f}"

    # 6. Lagna Point (good currency only, no debt)
    _la_lagna_sim = sim_good_total - sim_bad_total
    _la_lagna_pt_score = _la_lagna_sim
    _la_lagna_pt_notes = f"Sim Net = {_la_lagna_sim:.2f}"

    # 6b. Navamsa Lagna Score (standalone)
    _la_nav_score = sim_nav_net_score
    _la_nav_notes = f"Navamsa Lagna Net Score = {sim_nav_net_score:.2f}"

    # 7. Sun: (Maraivu adj Strength + Maraivu adj Score) / 2
    _la_sun_adj_str = _planet_maraivu_adj_strengths.get('Sun', 0.0)
    _la_sun_adj_nps = _nps_score_dict.get('Sun_adjusted', 0.0)
    _la_sun_suchama = _suchama_score_dict.get('Sun', 0.0)
    _la_sun_raw = (_la_sun_adj_str + _la_sun_adj_nps) / 2.0
    _la_sun_score = _la_sun_raw
    _la_sun_notes = f"(Str {_la_sun_adj_str:.2f} + AdjNPS {_la_sun_adj_nps:.2f})/2 = {_la_sun_raw:.2f}"

    # 8. 9th House Points
    _la_h9_raw = _house_total_points.get(9, 0.0)
    _la_h9_score = _la_h9_raw
    _la_h9_notes = f"House 9 total HP = {_la_h9_raw:.2f}"

    # 9. AG Bonus: weighted sum (LLStr includes Suchama)
    _ag_ll_str_combined = _la_ll_str_score + _la_ll_suchama_score
    _ag_moon   = _la_moon_score * 25.0 / 100.0
    _ag_ll     = _la_ll_score * 12.5 / 100.0
    _ag_ll_str = _ag_ll_str_combined * 12.5 / 100.0
    _ag_h1     = _la_h1_score * 40.0 / 100.0
    _ag_lp     = _la_lagna_pt_score * 10.0 / 100.0
    _ag_nav    = _la_nav_score * 5.0 / 100.0
    _ag_total  = _ag_moon + _ag_ll + _ag_ll_str + _ag_h1 + _ag_lp + _ag_nav
    _ag_notes  = (f"Moon({_la_moon_score:.2f}*25%)={_ag_moon:.2f} + "
                  f"LL({_la_ll_score:.2f}*12.5%)={_ag_ll:.2f} + "
                  f"LLStr+Suchama({_la_ll_str_score:.2f}+{_la_ll_suchama_score:.2f}={_ag_ll_str_combined:.2f}*12.5%)={_ag_ll_str:.2f} + "
                  f"H1({_la_h1_score:.2f}*40%)={_ag_h1:.2f} + "
                  f"LP({_la_lagna_pt_score:.2f}*10%)={_ag_lp:.2f} + "
                  f"NavLagna({_la_nav_score:.2f}*5%)={_ag_nav:.2f}")

    # 10. Bhuvi Bonus: weighted sum (LLStr includes Suchama)
    _bv_ll_str_combined = _la_ll_str_score + _la_ll_suchama_score
    _bv_moon   = _la_moon_score * 20.0 / 100.0
    _bv_ll     = _la_ll_score * 10.0 / 100.0
    _bv_ll_str = _bv_ll_str_combined * 10.0 / 100.0
    _bv_h1     = _la_h1_score * 30.0 / 100.0
    _bv_lp     = _la_lagna_pt_score * 5.0 / 100.0
    _bv_nav    = _la_nav_score * 5.0 / 100.0
    _bv_sun    = _la_sun_score * 10.0 / 100.0
    _bv_h9     = _la_h9_score * 10.0 / 100.0
    _bv_total  = _bv_moon + _bv_ll + _bv_ll_str + _bv_h1 + _bv_lp + _bv_nav + _bv_sun + _bv_h9
    _bv_notes  = (f"Moon({_la_moon_score:.2f}*20%)={_bv_moon:.2f} + "
                  f"LL({_la_ll_score:.2f}*10%)={_bv_ll:.2f} + "
                  f"LLStr+Suchama({_la_ll_str_score:.2f}+{_la_ll_suchama_score:.2f}={_bv_ll_str_combined:.2f}*10%)={_bv_ll_str:.2f} + "
                  f"H1({_la_h1_score:.2f}*30%)={_bv_h1:.2f} + "
                  f"LP({_la_lagna_pt_score:.2f}*5%)={_bv_lp:.2f} + "
                  f"NavLagna({_la_nav_score:.2f}*5%)={_bv_nav:.2f} + "
                  f"Sun({_la_sun_score:.2f}*10%)={_bv_sun:.2f} + "
                  f"H9({_la_h9_score:.2f}*10%)={_bv_h9:.2f}")

    lagna_analysis_rows = [
        ["Moon's Light",       f"{_la_moon_score:.2f}", _la_moon_notes],
        ['Lagna Lord Score',   f"{_la_ll_score:.2f}",   _la_ll_notes],
        ['Lagna Lord Strength',f"{_la_ll_str_score:.2f}", _la_ll_str_notes],
        ['Lagna Lord Suchama', f"{_la_ll_suchama_score:.2f}", _la_ll_suchama_notes],
        ['1st House Points',   f"{_la_h1_score:.2f}",   _la_h1_notes],
        ['Lagna Point',        f"{_la_lagna_pt_score:.2f}", _la_lagna_pt_notes],
        ['Navamsa Lagna Score',f"{_la_nav_score:.2f}",  _la_nav_notes],
        ['Sun Score',          f"{_la_sun_score:.2f}",  _la_sun_notes],
        ['9th House Points',   f"{_la_h9_score:.2f}",   _la_h9_notes],
        ['AG Bonus',           f"{_ag_total:.2f}",      _ag_notes],
        ['Bhuvi Bonus',        f"{_bv_total:.2f}",      _bv_notes],
    ]
    df_lagna_analysis = pd.DataFrame(lagna_analysis_rows, columns=['Metric', 'Score (out of 100)', 'Notes'])
    # ====== END LAGNA ANALYSIS TABLE ======

    return {
        'name': name, 'df_planets': df_planets, 'df_navamsa_exchange': df_navamsa_exchange,
        'df_navamsa_phase2': df_navamsa_phase2,
        'df_navamsa_phase3': df_navamsa_phase3,
        'df_phase1': df_phase1, 'df_phase2': df_phase2, 'df_phase3': df_phase3, 'df_phase4': df_phase4,
        'df_phase5': df_phase5, 'df_leftover_aspects': df_leftover_aspects,
        'df_jupiter_poison_notes': df_jupiter_poison_notes,
        'df_house_reserves': df_house_reserves,
        'df_bonus': df_bonus,
        'df_lagna_analysis': df_lagna_analysis,
        'df_normalized_planet_scores': df_normalized_planet_scores,
        'df_house_points': df_house_points,
        'df_planet_strengths': df_planet_strengths,
        'df_rasi': df_rasi, 'df_nav': df_nav,
        'df_house_status': df_house_status, 'dasa_periods_filtered': dasa_filtered,
        'lagna_sid': lagna_sid, 'nav_lagna': nav_lagna, 'lagna_sign': lagna_sign,
        'nav_lagna_sign': get_sign(nav_lagna), 'moon_rasi': get_sign(moon_lon),
        'moon_nakshatra': get_nakshatra_details(moon_lon)[0], 'moon_pada': get_nakshatra_details(moon_lon)[1],
        'selected_depth': depth_map[max_depth], 'utc_dt': utc_dt, 'max_depth': max_depth,
        'house_to_planets_rasi': house_planets_rasi, 'house_to_planets_nav': house_planets_nav
    }

# South Indian plotter
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

# Streamlit UI
st.set_page_config(page_title="Buvi Horoscope", layout="wide")
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

st.title("Buvi Astrology Data Generator")

if 'chart_data' not in st.session_state: st.session_state.chart_data = None
if 'search_results' not in st.session_state: st.session_state.search_results = []

@st.cache_resource
def get_geolocator():
    geolocator = Nominatim(user_agent="vedic_astro_app")
    return RateLimiter(geolocator.geocode, min_delay_seconds=1)

geocode = get_geolocator()
_tf = TimezoneFinder()

def tz_for_latlon(lat, lon):
    tzname = _tf.timezone_at(lng=lon, lat=lat)
    if not tzname: return pytz.UTC
    return pytz.timezone(tzname)

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

st.subheader("Birth Details")
name = st.text_input("Name", placeholder="Enter full name")
c1, c2 = st.columns(2)
with c1:
    birth_date = st.date_input("Birth Date", value=datetime.now().date(),
                               min_value=datetime(1,1,1).date(), max_value=datetime(2200,12,31).date())
with c2:
    birth_time = st.text_input("Birth Time (HH:MM in 24-hour format)", placeholder="14:30")

use_custom_coords = st.checkbox("Custom birth latitude and longitude?")
if use_custom_coords:
    clat, clon = st.columns(2)
    with clat: lat = st.number_input("Birth Latitude", value=13.08, format="%.4f")
    with clon: lon = st.number_input("Birth Longitude", value=80.27, format="%.4f")
else:
    birth_city_query = st.text_input("Birth City", value="Chennai", placeholder="Start typing birth city name...", key="birth_city_input")
    if birth_city_query and len(birth_city_query) >= 2:
        try:
            locations = geocode(birth_city_query, exactly_one=False, limit=5)
            st.session_state.search_results = [{'display': loc.address, 'lat': loc.latitude, 'lon': loc.longitude, 'address': loc.address} for loc in (locations or [])]
        except: st.session_state.search_results = []
    else: st.session_state.search_results = []

    if st.session_state.search_results:
        opts = [r['display'] for r in st.session_state.search_results]
        sel = st.selectbox("Select birth location", options=opts)
        i = opts.index(sel)
        lat = st.session_state.search_results[i]['lat']; lon = st.session_state.search_results[i]['lon']
    else:
        city_key = (birth_city_query or "").title()
        if city_key in cities_fallback:
            lat = cities_fallback[city_key]['lat']; lon = cities_fallback[city_key]['lon']
        else: lat, lon = 13.08, 80.27

# Auto-detect timezone from lat/lon and birth date
def _compute_tz_offset(lat_val, lon_val, date_obj):
    """Compute timezone UTC offset in hours for the given lat/lon and date."""
    try:
        tz = tz_for_latlon(lat_val, lon_val)
        # Use noon on the birth date to determine the UTC offset (handles DST correctly)
        naive_dt = datetime.combine(date_obj, datetime.min.time().replace(hour=12))
        localized_dt = tz.localize(naive_dt)
        offset_seconds = localized_dt.utcoffset().total_seconds()
        return offset_seconds / 3600.0, tz.zone
    except:
        return 5.5, "Asia/Kolkata"

auto_tz_offset, auto_tz_name = _compute_tz_offset(lat, lon, birth_date)

st.info(f"📍 Lat: {lat:.4f}, Lon: {lon:.4f} → Timezone: **{auto_tz_name}** (UTC {'+' if auto_tz_offset >= 0 else ''}{auto_tz_offset:g}h)")

override_tz = st.checkbox("Override auto-detected timezone?")
if override_tz:
    tz_offset = st.number_input("Timezone offset at birth (hrs)", value=auto_tz_offset, step=0.5)
else:
    tz_offset = auto_tz_offset

max_depth_options = {1:'Dasa only',2:'Dasa + Bhukti',3:'Dasa + Bhukti + Anthara',4:'Dasa + Bhukti + Anthara + Sukshma',5:'Dasa + Bhukti + Anthara + Sukshma + Prana',6:'Dasa + Bhukti + Anthara + Sukshma + Prana + Sub-Prana'}
selected_depth_str = st.selectbox("Generate up to (depth)", list(max_depth_options.values()), index=3)
max_depth = [k for k,v in max_depth_options.items() if v == selected_depth_str][0]

if st.button("Generate Chart", use_container_width=True):
    if not name or not birth_time: st.error("Please enter Name and Birth Time.")
    else:
        try:
            with st.spinner("Calculating chart..."):
                st.session_state.chart_data = compute_chart(name, birth_date, birth_time, lat, lon, tz_offset, max_depth)
            st.rerun()
        except Exception as e: st.error(f"Error: {e}")

def _serialize_dasa_periods(periods, level_names=None, depth=0):
    """Recursively serialize Dasa/Bhukti/Antara periods to JSON-friendly dicts."""
    if level_names is None:
        level_names = ['Dasa', 'Bhukti', 'Anthara', 'Sukshma', 'Prana', 'Sub-Prana']
    result = []
    current_level = level_names[depth] if depth < len(level_names) else f"Level-{depth+1}"
    for lord, start, end, subs in periods:
        entry = {
            'level': current_level,
            'planet': lord,
            'start': start.strftime('%Y-%m-%d %H:%M'),
            'end': end.strftime('%Y-%m-%d %H:%M'),
            'duration': duration_str(end - start, current_level.lower())
        }
        if subs:
            entry['sub_periods'] = _serialize_dasa_periods(subs, level_names, depth + 1)
        result.append(entry)
    return result

def build_export_json(cd):
    """Build a comprehensive JSON dict from chart data, excluding Phase 1-4
    currency exchange and Navamsa Phase 1-3."""
    data = {}

    # ── 1. Chart Summary ──
    data['chart_summary'] = {
        'name': cd['name'],
        'lagna_sign': cd['lagna_sign'],
        'lagna_degrees': round(cd['lagna_sid'], 2),
        'navamsa_lagna_sign': cd['nav_lagna_sign'],
        'navamsa_lagna_degrees': round(cd['nav_lagna'], 2),
        'moon_rasi': cd['moon_rasi'],
        'moon_nakshatra': cd['moon_nakshatra'],
        'moon_pada': cd['moon_pada'],
        'dasa_depth': cd['selected_depth'],
        'utc_datetime': cd['utc_dt'].strftime('%Y-%m-%d %H:%M:%S')
    }

    # ── 2. Planetary Positions ──
    data['planetary_positions'] = cd['df_planets'].to_dict(orient='records')

    # ── 3. Rasi Chart (D1) ──
    data['rasi_chart'] = cd['df_rasi'].to_dict(orient='records')

    # ── 4. Navamsa Chart (D9) ──
    data['navamsa_chart'] = cd['df_nav'].to_dict(orient='records')

    # ── 5. House Analysis ──
    data['house_analysis'] = cd['df_house_status'].to_dict(orient='records')

    # ── 6. House Bonus Points (Reserve) ──
    data['house_bonus_reserves'] = cd['df_house_reserves'].to_dict(orient='records')

    # ── 7. Currency Exchange Phase 5 (Final) ──
    data['currency_exchange_phase5'] = cd['df_phase5'].to_dict(orient='records')

    # ── 8. Leftover Aspect Clones (Phase 5) ──
    data['leftover_aspect_clones'] = cd['df_leftover_aspects'].to_dict(orient='records')

    # ── 9. Jupiter Poison Diagnostic ──
    data['jupiter_poison_diagnostic'] = cd['df_jupiter_poison_notes'].to_dict(orient='records')

    # ── 10. Normalized Planet Scores ──
    data['normalized_planet_scores'] = cd['df_normalized_planet_scores'].to_dict(orient='records')

    # ── 11. House Points Analysis ──
    data['house_points_analysis'] = cd['df_house_points'].to_dict(orient='records')

    # ── 12. Planet Strengths ──
    data['planet_strengths'] = cd['df_planet_strengths'].to_dict(orient='records')

    # ── 13. Lagna Point Score Simulation ──
    data['lagna_simulation'] = cd['df_bonus'].to_dict(orient='records')

    # ── 14. Lagna Analysis ──
    data['lagna_analysis'] = cd['df_lagna_analysis'].to_dict(orient='records')

    # ── 15. Vimshottari Dasa (Full Hierarchy: Dasa → Bhukti → Anthara ...) ──
    data['vimshottari_dasa'] = _serialize_dasa_periods(cd['dasa_periods_filtered'])

    return data

def show_png(fig):
    fig.tight_layout(pad=0.10)
    st.pyplot(fig, use_container_width=False, dpi=300)

if st.session_state.chart_data:
    cd = st.session_state.chart_data
    st.markdown("---")
    st.markdown(f"""
    <div class="summary-box">
        <h3>Chart Summary</h3>
        <div class="summary-item"><strong>Name:</strong> {cd['name']}</div>
        <div class="summary-item"><strong>Lagna:</strong> {cd['lagna_sign']} ({cd['lagna_sid']:.2f}deg)</div>
        <div class="summary-item"><strong>Rasi (Moon Sign):</strong> {cd['moon_rasi']}</div>
        <div class="summary-item"><strong>Nakshatra:</strong> {cd['moon_nakshatra']} (Pada {cd['moon_pada']})</div>
    </div>
    """, unsafe_allow_html=True)

    st.subheader("Planetary Positions")
    st.dataframe(cd['df_planets'], hide_index=True, use_container_width=True)

    st.subheader("House Bonus Points (Reserve)")
    st.dataframe(cd['df_house_reserves'], hide_index=True, use_container_width=True)

    st.subheader("Navamsa Exchange (Phase 1)")
    st.dataframe(cd['df_navamsa_exchange'], hide_index=True, use_container_width=True)

    st.subheader("Navamsa Phase 2 Exchange")
    st.dataframe(cd['df_navamsa_phase2'], hide_index=True, use_container_width=True)

    st.subheader("Navamsa Phase 3 Exchange")
    st.dataframe(cd['df_navamsa_phase3'], hide_index=True, use_container_width=True)

    st.subheader("Currency Exchange Phase 1")
    st.dataframe(cd['df_phase1'], hide_index=True, use_container_width=True)

    st.subheader("Currency Exchange Phase 2")
    st.dataframe(cd['df_phase2'], hide_index=True, use_container_width=True)

    st.subheader("Currency Exchange Phase 3")
    st.dataframe(cd['df_phase3'], hide_index=True, use_container_width=True)

    st.subheader("Currency Exchange Phase 4")
    st.dataframe(cd['df_phase4'], hide_index=True, use_container_width=True)

    st.subheader("Currency Exchange Phase 5")
    st.dataframe(cd['df_phase5'], hide_index=True, use_container_width=True)

    st.subheader("Leftover Aspect Clones (Phase 5)")
    st.dataframe(cd['df_leftover_aspects'], hide_index=True, use_container_width=True)

    st.subheader("Jupiter Poison Diagnostic")
    st.dataframe(cd['df_jupiter_poison_notes'], hide_index=True, use_container_width=True)

    st.subheader("Normalized Planet Scores")
    st.dataframe(cd['df_normalized_planet_scores'], hide_index=True, use_container_width=True)

    st.subheader("House Points Analysis")
    st.dataframe(cd['df_house_points'], hide_index=True, use_container_width=True)

    st.subheader("Planet Strengths")
    st.dataframe(cd['df_planet_strengths'], hide_index=True, use_container_width=True)

    st.subheader("Lagna Point Score Simulation")
    st.dataframe(cd['df_bonus'], hide_index=True, use_container_width=True)

    st.subheader("Lagna Analysis")
    st.dataframe(cd['df_lagna_analysis'], hide_index=True, use_container_width=True)

    st.subheader("Rasi (D1) & Navamsa (D9) - South Indian")
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
    dasa_rows = [{'Planet': lord, 'Start': s.strftime('%Y-%m-%d'), 'End': e.strftime('%Y-%m-%d'), 'Duration': duration_str(e-s,'dasa')} for lord, s, e, _ in cd['dasa_periods_filtered']]
    st.dataframe(pd.DataFrame(dasa_rows), hide_index=True, use_container_width=True)

    dp = cd['dasa_periods_filtered']
    if cd['max_depth'] >= 2:
        with st.expander("View Sub-periods (Bhukti / Anthara / Sukshma)", expanded=False):
            # --- Bhukti level ---
            d_opt = [f"{p[0]} ({p[1].strftime('%Y-%m-%d')} → {p[2].strftime('%Y-%m-%d')})" for p in dp]
            sel_dasa = st.selectbox("Select Dasa:", d_opt, key="sel_dasa")
            sel_dasa_idx = d_opt.index(sel_dasa)
            bhuktis = dp[sel_dasa_idx][3]
            if bhuktis:
                st.markdown("**Bhukti (Sub-periods)**")
                st.dataframe(pd.DataFrame([{'Planet': l, 'Start': s.strftime('%Y-%m-%d'), 'End': e.strftime('%Y-%m-%d'), 'Duration': duration_str(e-s,'bhukti')} for l,s,e,_ in bhuktis]), hide_index=True, use_container_width=True)

            # --- Anthara level ---
            if cd['max_depth'] >= 3 and bhuktis:
                b_opt = [f"{p[0]} ({p[1].strftime('%Y-%m-%d')} → {p[2].strftime('%Y-%m-%d')})" for p in bhuktis]
                sel_bhukti = st.selectbox("Select Bhukti to view Anthara:", b_opt, key="sel_bhukti")
                sel_bhukti_idx = b_opt.index(sel_bhukti)
                antaras = bhuktis[sel_bhukti_idx][3]
                if antaras:
                    st.markdown("**Anthara (Sub-sub-periods)**")
                    st.dataframe(pd.DataFrame([{'Planet': l, 'Start': s.strftime('%Y-%m-%d'), 'End': e.strftime('%Y-%m-%d'), 'Duration': duration_str(e-s,'anthara')} for l,s,e,_ in antaras]), hide_index=True, use_container_width=True)

                # --- Sukshma level ---
                if cd['max_depth'] >= 4 and antaras:
                    a_opt = [f"{p[0]} ({p[1].strftime('%Y-%m-%d')} → {p[2].strftime('%Y-%m-%d')})" for p in antaras]
                    sel_antara = st.selectbox("Select Anthara to view Sukshma:", a_opt, key="sel_antara")
                    sel_antara_idx = a_opt.index(sel_antara)
                    sukshmas = antaras[sel_antara_idx][3]
                    if sukshmas:
                        st.markdown("**Sukshma (Sub-sub-sub-periods)**")
                        st.dataframe(pd.DataFrame([{'Planet': l, 'Start': s.strftime('%Y-%m-%d %H:%M'), 'End': e.strftime('%Y-%m-%d %H:%M'), 'Duration': duration_str(e-s,'sukshma')} for l,s,e,_ in sukshmas]), hide_index=True, use_container_width=True)

    st.subheader("Current City - Live Micro-Periods")
    current_city_query = st.text_input("Enter your CURRENT city", placeholder="e.g., Chennai", key="current_city_input")
    depth_choice = st.selectbox("Depth to inspect", ["Sukshma", "Prana", "Sub-Prana"])
    
    if st.button("Show current micro-periods", use_container_width=True):
        if current_city_query:
            try:
                cur_locs = geocode(current_city_query, exactly_one=False, limit=1)
                if cur_locs:
                    cur = cur_locs[0]; tz = tz_for_latlon(cur.latitude, cur.longitude)
                    now_local = datetime.now(tz); now_utc_naive = now_local.astimezone(pytz.UTC).replace(tzinfo=None)
                    active_path = find_active_path_to_depth(dp, now_utc_naive, _DEPTH_NAME_TO_INT[depth_choice])
                    flat_at_depth = collect_periods_at_depth(dp, _DEPTH_NAME_TO_INT[depth_choice])
                    
                    st.success(f"Time zone: {tz.zone} | Local now: {now_local.strftime('%Y-%m-%d %H:%M')}")
                    if active_path:
                        tbl = []
                        idx_found = -1
                        for i, (lord, s, e) in enumerate(flat_at_depth):
                            if s <= now_utc_naive < e:
                                idx_found = i
                                break
                        if idx_found != -1:
                            for l,st_t,en_t in flat_at_depth[idx_found : idx_found+6]:
                                tbl.append({"Lord": l, "Start (local)": st_t.replace(tzinfo=pytz.UTC).astimezone(tz).strftime('%Y-%m-%d %H:%M'), "End (local)": en_t.replace(tzinfo=pytz.UTC).astimezone(tz).strftime('%Y-%m-%d %H:%M'), "Duration": duration_str(en_t-st_t, depth_choice.lower())})
                        st.dataframe(pd.DataFrame(tbl), hide_index=True, use_container_width=True)
            except Exception as e: st.error(f"Error: {e}")

    # ====== EXPORT ALL DATA AS JSON ======
    st.markdown("---")
    st.subheader("Export All Data as JSON")
    st.caption("Excludes: Phase 1-4 Currency Exchange & Navamsa Phase 1-3. Includes everything else + full Dasa/Bhukti/Antara hierarchy.")
    if st.button("Generate JSON", use_container_width=True, key="gen_json_btn"):
        export_data = build_export_json(cd)
        json_str = json.dumps(export_data, indent=2, ensure_ascii=False)
        st.session_state['export_json'] = json_str

    if 'export_json' in st.session_state and st.session_state['export_json']:
        json_str = st.session_state['export_json']

        st.download_button(
            label="Download JSON File",
            data=json_str,
            file_name=f"{cd['name'].replace(' ', '_')}_chart_data.json",
            mime="application/json",
            use_container_width=True
        )

        with st.expander("View / Copy JSON", expanded=False):
            st.code(json_str, language="json")

else: st.info("Enter birth details above and click 'Generate Chart' to begin")

st.markdown("---")
st.caption("Buvi Astrology Data Generator")
