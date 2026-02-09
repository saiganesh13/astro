import streamlit as st
from datetime import datetime, timedelta
from math import sin, cos, tan, atan2, degrees, radians
from collections import defaultdict
import pandas as pd
from geopy.geocoders import Nominatim
from geopy.extra.rate_limiter import RateLimiter
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from timezonefinder import TimezoneFinder
import pytz
import copy

# ---- Swiss Ephemeris (primary) with Astropy fallback ----
try:
    import swisseph as swe
    swe.set_sid_mode(swe.SIDM_LAHIRI)
    USE_SWISSEPH = True
except ImportError:
    USE_SWISSEPH = False
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

# ---- Astropy fallback functions (only used when swisseph unavailable) ----
def get_lahiri_ayanamsa_fallback(year):
    """Improved Lahiri ayanamsa - polynomial fit accurate across centuries."""
    T = (year - 2000) / 100.0
    # More accurate polynomial coefficients for Lahiri ayanamsa
    ayan = 23.853 + 1.39722 * T + 0.000306 * T**2 + 0.000016 * T**3
    return ayan % 360

def get_obliquity(d):
    T = d/36525.0
    return ((((-4.34e-8*T - 5.76e-7)*T + 0.0020034)*T - 1.831e-4)*T - 46.836769)*T/3600 + 23.4392794444444

def get_gmst(d):
    T = d/36525.0
    return (67310.54841 + (3155760000 + 8640184.812866)*T + 0.093104*T**2 - 6.2e-6*T**3)/3600 % 24

def get_ascendant_fallback(jd, lat, lon):
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
    if 'Good' in c_key:
        return True
    return False

def is_sun_or_moon_currency(c_key):
    return 'Sun' in c_key or 'Moon' in c_key

# ---- Swiss Ephemeris planet computation ----
def compute_positions_swisseph(utc_dt, lat, lon):
    """Compute all sidereal planet positions + ascendant using Swiss Ephemeris."""
    year = utc_dt.year
    month = utc_dt.month
    day = utc_dt.day
    hour = utc_dt.hour + utc_dt.minute / 60.0 + utc_dt.second / 3600.0
    
    # Determine calendar type: Julian before Oct 15, 1582; Gregorian after
    if utc_dt < datetime(1582, 10, 15):
        cal_type = swe.JUL_CAL
    else:
        cal_type = swe.GREG_CAL
    
    jd = swe.julday(year, month, day, hour, cal_type)
    
    # Set Lahiri ayanamsa
    swe.set_sid_mode(swe.SIDM_LAHIRI)
    ayan = swe.get_ayanamsa_ut(jd)
    
    # Planet map: name -> swisseph ID
    planet_ids = {
        'sun': swe.SUN, 'moon': swe.MOON, 'mercury': swe.MERCURY,
        'venus': swe.VENUS, 'mars': swe.MARS, 'jupiter': swe.JUPITER,
        'saturn': swe.SATURN
    }
    
    # Compute sidereal longitudes
    lon_sid = {}
    for name, pid in planet_ids.items():
        result = swe.calc_ut(jd, pid, swe.FLG_SIDEREAL)
        lon_sid[name] = result[0][0] % 360
    
    # Rahu (Mean Node) - tropical, then convert to sidereal
    rahu_trop = swe.calc_ut(jd, swe.MEAN_NODE)[0][0]
    lon_sid['rahu'] = (rahu_trop - ayan) % 360
    lon_sid['ketu'] = (lon_sid['rahu'] + 180) % 360
    
    # Ascendant - tropical from swe.houses, then convert to sidereal
    cusps, ascmc = swe.houses(jd, lat, lon, b'P')  # Placidus
    asc_trop = ascmc[0]
    lagna_sid = (asc_trop - ayan) % 360
    
    return lon_sid, lagna_sid, ayan, jd

# ---- Astropy fallback planet computation ----
def compute_positions_astropy(utc_dt, lat, lon):
    """Fallback: compute positions using astropy."""
    from astropy.time import Time
    from astropy.coordinates import get_body, solar_system_ephemeris, GeocentricTrueEcliptic
    
    t = Time(utc_dt)
    jd = t.jd
    ayan = get_lahiri_ayanamsa_fallback(utc_dt.year)
    
    with solar_system_ephemeris.set('builtin'):
        lon_trop = {}
        for nm in ['sun','moon','mercury','venus','mars','jupiter','saturn']:
            ecl = get_body(nm, t).transform_to(GeocentricTrueEcliptic())
            lon_trop[nm] = ecl.lon.deg
    
    d = jd - 2451545.0; T = d/36525.0
    omega = (125.04452 - 1934.136261*T + 0.0020708*T**2 + T**3/450000) % 360
    lon_trop['rahu'] = omega
    lon_trop['ketu'] = (omega + 180) % 360
    
    lon_sid = {p: get_sidereal_lon(lon_trop[p], ayan) for p in lon_trop}
    lagna_sid = get_sidereal_lon(get_ascendant_fallback(jd, lat, lon), ayan)
    
    return lon_sid, lagna_sid, ayan, jd


def compute_chart(name, date_obj, time_str, lat, lon, tz_offset, max_depth):
    try:
        hour, minute = map(int, time_str.split(':'))
        if not (0<=hour<=23 and 0<=minute<=59): raise ValueError
    except:
        raise ValueError("Time must be in HH:MM format (24-hour)")
    
    local_dt = datetime.combine(date_obj, datetime.min.time().replace(hour=hour, minute=minute))
    utc_dt = local_dt - timedelta(hours=tz_offset)
    
    # ---- Compute planet positions using best available engine ----
    if USE_SWISSEPH:
        lon_sid, lagna_sid, ayan, jd = compute_positions_swisseph(utc_dt, lat, lon)
    else:
        lon_sid, lagna_sid, ayan, jd = compute_positions_astropy(utc_dt, lat, lon)
    
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
    
    # Navamsa Exchange Cycle
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
                
                if debtor == 'Ketu' and t_name not in ['Sun', 'Moon']: continue
                
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
        p3_cycle_limit = 200
        p3_cycles = 0
        
        p3_standard_malefics = ['Saturn', 'Rahu', 'Ketu', 'Mars', 'Sun']
        p3_standard_benefics = ['Jupiter', 'Venus', 'Mercury']
        p3_malefic_debtor_order = ['Rahu', 'Sun', 'Saturn', 'Mars', 'Ketu']
        
        while p3_cycles < p3_cycle_limit:
            p3_cycles += 1
            p3_something_happened = False
            
            if house_11_pot <= 0.001:
                break
            
            moon_bad_currency_p3 = phase3_data['Moon']['p3_inventory'].get('Bad Moon', 0.0)
            moon_is_malefic_p3 = moon_bad_currency_p3 > 0.001
            
            house_11_malefics = []
            for p in planets_in_house_11:
                if p in p3_standard_malefics:
                    house_11_malefics.append(p)
                elif p == 'Moon' and moon_is_malefic_p3:
                    house_11_malefics.append(p)
            
            house_11_benefics = []
            for p in planets_in_house_11:
                if p in p3_standard_benefics:
                    house_11_benefics.append(p)
                elif p == 'Moon' and not moon_is_malefic_p3:
                    house_11_benefics.append(p)
            
            def get_p3_malefic_rank(p):
                if p == 'Moon':
                    vol = phase3_data['Moon']['volume']
                    if vol > 0:
                        bad_pct = (moon_bad_currency_p3 / vol) * 100
                    else:
                        bad_pct = 0
                    if bad_pct > 25:
                        return 2.5
                    else:
                        return 3.5
                else:
                    if p in p3_malefic_debtor_order:
                        return p3_malefic_debtor_order.index(p)
                    return 99
            
            house_11_malefics_sorted = sorted(house_11_malefics, key=get_p3_malefic_rank)
            
            def get_p3_benefic_debt_pct(p):
                vol = phase3_data[p]['volume']
                debt = abs(phase3_data[p]['p3_current_debt'])
                if vol > 0:
                    return (debt / vol) * 100
                return 0
            
            house_11_benefics_sorted = sorted(house_11_benefics, key=lambda p: -get_p3_benefic_debt_pct(p))
            
            for malefic in house_11_malefics_sorted:
                if house_11_pot <= 0.001:
                    break
                
                debt = phase3_data[malefic]['p3_current_debt']
                
                if debt < -0.001:
                    needed = abs(debt)
                    take = min(1.0, needed, house_11_pot)
                    
                    if take > 0.001:
                        house_11_pot -= take
                        phase3_data[malefic]['p3_inventory']['Good Moon'] += take
                        phase3_data[malefic]['p3_current_debt'] += take
                        p3_something_happened = True
            
            all_malefics_cleared = True
            for malefic in house_11_malefics:
                if phase3_data[malefic]['p3_current_debt'] < -0.001:
                    all_malefics_cleared = False
                    break
            
            if all_malefics_cleared or len(house_11_malefics) == 0:
                for benefic in house_11_benefics_sorted:
                    if house_11_pot <= 0.001:
                        break
                    
                    debt = phase3_data[benefic]['p3_current_debt']
                    
                    if debt < -0.001:
                        needed = abs(debt)
                        take = min(1.0, needed, house_11_pot)
                        
                        if take > 0.001:
                            house_11_pot -= take
                            phase3_data[benefic]['p3_inventory']['Good Moon'] += take
                            phase3_data[benefic]['p3_current_debt'] += take
                            p3_something_happened = True
            
            if not p3_something_happened:
                break
    
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
    
    for sign_name_key, (gifter, multiplier) in gift_pot_config.items():
        gifter_sthana = planet_data[gifter]['sthana']
        pot_value = multiplier * (gifter_sthana / 100.0)
        pot_inventory[sign_name_key] = pot_value
        pot_currency_type[sign_name_key] = gifter
    
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
    
    def get_p5_currency_rank_score(c_key):
        if c_key == 'Jupiter': return 990
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
    
    for current_planet in PLANET_SEQUENCE:
        if current_planet not in ASPECT_RULES:
            continue
        
        aspect_offsets = ASPECT_RULES[current_planet]
        parent_data = phase5_data[current_planet]
        parent_L = parent_data['L']
        parent_debt = parent_data['p5_current_debt']
        
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
                
                clone_value = aspect_pct * good_sum
                if clone_value > 0.001:
                    clone_inventory['Good Saturn'] = clone_value
                
                clone_debt = parent_debt * aspect_pct
                clone_type = 'Active'
                
            elif current_planet == 'Mars':
                good_sum = 0.0
                for k, v in parent_inv.items():
                    if is_good_currency(k) and v > 0.001:
                        if k == 'Good Mars':
                            good_sum += v
                        else:
                            good_sum += v / 2.0
                
                clone_value = aspect_pct * good_sum
                if clone_value > 0.001:
                    clone_inventory['Good Mars'] = clone_value
                
                clone_debt = parent_debt * aspect_pct
                clone_type = 'Active'
                
            elif current_planet in ['Jupiter', 'Venus', 'Mercury']:
                good_sum = 0.0
                for k, v in parent_inv.items():
                    if is_good_currency(k) and v > 0.001:
                        good_sum += v
                
                clone_value = aspect_pct * good_sum
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
                clone_value = aspect_pct * total_value
                if clone_value > 0.001:
                    clone_inventory['Good Moon'] = clone_value
                
                clone_debt = 0.0
                clone_type = 'Passive'
            
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
        
        # Part B: The Interaction Cycle
        p5_cycle_limit = 500
        p5_cycles = 0
        
        while p5_cycles < p5_cycle_limit:
            p5_cycles += 1
            p5_something_happened = False
            
            for clone in clones:
                clone_L = clone['L']
                
                # Step 1 - Virtual Malefic Clones Pull (Active Pulling)
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
                
                # Step 2 - Real Malefics Pull (from Clone's original inventory only)
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
                inv_str,
                debt_str
            ])
    
    # Format Phase 5 Output - NOW WITH GOOD SCORE COLUMN
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
        
        # GOOD SCORE: 100 if no debt, otherwise 100 - abs(debt), minimum 0
        if d_val >= -0.01:
            phase5_data[p]['good_score'] = 100.0
        else:
            phase5_data[p]['good_score'] = max(0.0, 100.0 + d_val)
    
    for p in ['Sun','Moon','Mars','Mercury','Jupiter','Venus','Saturn','Rahu','Ketu']:
        d_p5 = phase5_data[p]
        phase5_rows.append([p, d_p5['currency_p5'], d_p5['debt_p5'], f"{d_p5['good_score']:.1f}"])
    
    df_phase5 = pd.DataFrame(phase5_rows, columns=['Planet', 'Currency [Phase 5]', 'Debt [Phase 5]', 'Good Score'])
    
    df_leftover_aspects = pd.DataFrame(leftover_aspects, columns=['Source Planet', 'Aspect Angle', 'Remaining Inventory', 'Final Debt'])
    
    # CREATE HOUSE RESERVES DATAFRAME
    reserve_rows = []
    for sign_name_iter in sign_names:
        currency_dict = house_reserves[sign_name_iter]
        
        currency_parts = []
        for currency_type, amount in currency_dict.items():
            if amount > 0.001:
                currency_parts.append(f"{currency_type}[{amount:.2f}]")
        
        reserve_str = ", ".join(currency_parts) if currency_parts else "-"
        reserve_rows.append([sign_name_iter, reserve_str])
    
    df_house_reserves = pd.DataFrame(reserve_rows, columns=['House Sign', 'Unutilized Bonus Points'])
    
    df_planets = pd.DataFrame(rows, columns=['Planet','Deg','Sign','Nakshatra','Pada','Ld/SL','Vargothuva',
                                             'Parivardhana',
                                             'Dig Bala (%)','Sthana Bala (%)','Status','Updated Status',
                                             'Volume', 'Default Currencies', 'Debt'])

    df_rasi = pd.DataFrame([[f"House {h}", get_sign((lagna_sid+(h-1)*30)%360), 
                             ', '.join(sorted(house_planets_rasi[h])) if house_planets_rasi[h] else 'Empty'] 
                            for h in range(1,13)], columns=['House','Sign','Planets'])

    house_planets_nav = defaultdict(list)
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

    engine_used = "Swiss Ephemeris" if USE_SWISSEPH else "Astropy (fallback)"

    return {
        'name': name, 'df_planets': df_planets, 'df_navamsa_exchange': df_navamsa_exchange,
        'df_navamsa_phase2': df_navamsa_phase2,
        'df_navamsa_phase3': df_navamsa_phase3,
        'df_phase1': df_phase1, 'df_phase2': df_phase2, 'df_phase3': df_phase3, 'df_phase4': df_phase4,
        'df_phase5': df_phase5, 'df_leftover_aspects': df_leftover_aspects,
        'df_house_reserves': df_house_reserves,
        'df_rasi': df_rasi, 'df_nav': df_nav,
        'df_house_status': df_house_status, 'dasa_periods_filtered': dasa_filtered,
        'lagna_sid': lagna_sid, 'nav_lagna': nav_lagna, 'lagna_sign': lagna_sign,
        'nav_lagna_sign': get_sign(nav_lagna), 'moon_rasi': get_sign(moon_lon),
        'moon_nakshatra': get_nakshatra_details(moon_lon)[0], 'moon_pada': get_nakshatra_details(moon_lon)[1],
        'selected_depth': depth_map[max_depth], 'utc_dt': utc_dt, 'max_depth': max_depth,
        'house_to_planets_rasi': house_planets_rasi, 'house_to_planets_nav': house_planets_nav,
        'engine_used': engine_used
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

engine_label = "Swiss Ephemeris (high accuracy)" if USE_SWISSEPH else "Astropy (fallback)"
st.caption(f"Ephemeris engine: **{engine_label}** | Supported date range: **1 AD - 2200 AD**")

name = st.text_input("Name", placeholder="Enter full name")
c1, c2, c3 = st.columns(3)
with c1:
    # Extended date range: 1 AD to 2200 AD
    birth_date = st.date_input("Birth Date", value=datetime.now().date(),
                               min_value=datetime(1, 1, 1).date(), max_value=datetime(2200, 12, 31).date())
with c2:
    birth_time = st.text_input("Birth Time (HH:MM in 24-hour format)", placeholder="14:30")
with c3:
    tz_offset = st.number_input("Timezone offset at birth (hrs)", value=5.5, step=0.5)

use_custom_coords = st.checkbox("Custom birth latitude and longitude?")
if use_custom_coords:
    clat, clon = st.columns(2)
    with clat: lat = st.number_input("Birth Latitude", value=13.08, format="%.4f")
    with clon: lon = st.number_input("Birth Longitude", value=80.27, format="%.4f")
else:
    birth_city_query = st.text_input("Birth City", placeholder="Start typing birth city name...", key="birth_city_input")
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

max_depth_options = {1:'Dasa only',2:'Dasa + Bhukti',3:'Dasa + Bhukti + Anthara',4:'Dasa + Bhukti + Anthara + Sukshma',5:'Dasa + Bhukti + Anthara + Sukshma + Prana',6:'Dasa + Bhukti + Anthara + Sukshma + Prana + Sub-Prana'}
selected_depth_str = st.selectbox("Generate up to (depth)", list(max_depth_options.values()), index=2)
max_depth = [k for k,v in max_depth_options.items() if v == selected_depth_str][0]

if st.button("Generate Chart", use_container_width=True):
    if not name or not birth_time: st.error("Please enter Name and Birth Time.")
    else:
        try:
            with st.spinner("Calculating chart..."):
                st.session_state.chart_data = compute_chart(name, birth_date, birth_time, lat, lon, tz_offset, max_depth)
            st.rerun()
        except Exception as e: st.error(f"Error: {e}")

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
        <div class="summary-item"><strong>Calculation Engine:</strong> {cd['engine_used']}</div>
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

    st.subheader("Currency Exchange Phase 5 (Final)")
    st.dataframe(cd['df_phase5'], hide_index=True, use_container_width=True)

    st.subheader("Leftover Aspect Clones (Phase 5)")
    st.dataframe(cd['df_leftover_aspects'], hide_index=True, use_container_width=True)

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
        with st.expander("View Sub-periods", expanded=False):
            d_opt = [f"{p[0]} ({p[1].strftime('%Y-%m-%d')} - {p[2].strftime('%Y-%m-%d')})" for p in dp]
            sel = st.selectbox("Select Dasa:", d_opt)
            bhuktis = dp[d_opt.index(sel)][3]
            st.dataframe(pd.DataFrame([{'Planet': l, 'Start': s.strftime('%Y-%m-%d'), 'End': e.strftime('%Y-%m-%d'), 'Duration': duration_str(e-s,'bhukti')} for l,s,e,_ in bhuktis]), hide_index=True, use_container_width=True)

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

else: st.info("Enter birth details above and click 'Generate Chart' to begin")

st.markdown("---")
st.caption("Buvi Astrology Data Generator")
