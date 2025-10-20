import streamlit as st
from datetime import datetime, timedelta
from math import sin, cos, tan, atan2, degrees, radians
from astropy.time import Time
from astropy.coordinates import get_body, solar_system_ephemeris, GeocentricTrueEcliptic
from collections import defaultdict
import pandas as pd
from geopy.geocoders import Nominatim
from geopy.extra.rate_limiter import RateLimiter

# Note: Install required packages: pip install streamlit astropy geopy pandas

# City data fallback (Indian cities)
cities_fallback = {
    'Chennai': {'lat': 13.08, 'lon': 80.27},
    'Mumbai': {'lat': 19.07, 'lon': 72.88},
    'Delhi': {'lat': 28.61, 'lon': 77.23},
    'Bangalore': {'lat': 12.97, 'lon': 77.59},
    'Kolkata': {'lat': 22.57, 'lon': 88.36},
    'Hyderabad': {'lat': 17.39, 'lon': 78.49},
}

sign_names = ['Aries', 'Taurus', 'Gemini', 'Cancer', 'Leo', 'Virgo', 'Libra', 'Scorpio', 'Sagittarius', 'Capricorn', 'Aquarius', 'Pisces']

lords_full = ['Ketu', 'Venus', 'Sun', 'Moon', 'Mars', 'Rahu', 'Jupiter', 'Saturn', 'Mercury']
lords_short = ['Ke', 'Ve', 'Su', 'Mo', 'Ma', 'Ra', 'Ju', 'Sa', 'Me']
lords = lords_full

nak_names = [
    'Ashwini', 'Bharani', 'Krittika', 'Rohini', 'Mrigashira', 'Ardra', 'Punarvasu', 'Pushya', 'Ashlesha',
    'Magha', 'Purva Phalguni', 'Uttara Phalguni', 'Hasta', 'Chitra', 'Swati', 'Vishakha', 'Anuradha',
    'Jyeshta', 'Mula', 'Purva Ashadha', 'Uttara Ashadha', 'Shravana', 'Dhanishta', 'Shatabhisha',
    'Purva Bhadrapada', 'Uttara Bhadrapada', 'Revati'
]

years = [7, 20, 6, 10, 7, 18, 16, 19, 17] * 3

sign_lords = ['Mars', 'Venus', 'Mercury', 'Moon', 'Sun', 'Mercury', 'Venus', 'Mars', 'Jupiter', 'Saturn', 'Saturn', 'Jupiter']

def get_lahiri_ayanamsa(year):
    base = 23.853
    rate = 50.2388 / 3600.0
    return (base + (year - 2000) * rate) % 360

def get_obliquity(d):
    T = d / 36525.0
    oe = ((((-4.34e-8 * T - 5.76e-7) * T + 0.0020034) * T - 1.831e-4) * T - 46.836769) * T / 3600 + 23.4392794444444
    return oe

def get_gmst(d):
    T = d / 36525.0
    gmst = (67310.54841 + (3155760000 + 8640184.812866) * T + 0.093104 * T**2 - 6.2e-6 * T**3) / 3600 % 24
    return gmst

def get_ascendant(jd, lat, lon):
    d = jd - 2451545.0
    oe = get_obliquity(d)
    oer = radians(oe)
    gmst = get_gmst(d)
    lst = (gmst + lon / 15.0) % 24
    lstr = radians(lst * 15.0)
    sin_asc = cos(lstr)
    cos_asc = -(sin(lstr) * cos(oer) + tan(radians(lat)) * sin(oer))
    ascr = atan2(sin_asc, cos_asc)
    asc = degrees(ascr) % 360
    return asc

def get_sidereal_lon(trop_lon, ayan):
    return (trop_lon - ayan) % 360

def get_sign(lon):
    return sign_names[int(lon / 30)]

def get_house(lon, lagna_lon):
    s = int(lon / 30)
    l = int(lagna_lon / 30)
    return (s - l) % 12 + 1

def get_nakshatra_details(lon):
    deg_per_nak = 360 / 27
    nak_idx = int(lon // deg_per_nak) % 27
    nak_name = nak_names[nak_idx]
    pos_in_nak = lon % deg_per_nak
    deg_per_pada = deg_per_nak / 4
    pada = int(pos_in_nak // deg_per_pada) + 1
    star_lord_idx = nak_idx % 9
    star_lord = lords_short[star_lord_idx]
    fraction = pos_in_nak / deg_per_nak
    sub_num = int(fraction * 9)
    sub_lord_idx = (star_lord_idx + sub_num) % 9
    sub_lord = lords_short[sub_lord_idx]
    return nak_name, pada, star_lord, sub_lord

def generate_vimshottari_dasa(moon_lon):
    nak = int(moon_lon * 27 / 360)
    lord_idx = nak % 9
    y = years[lord_idx]
    full_nak_deg = 360 / 27
    pos_in_nak = moon_lon % full_nak_deg
    fraction = pos_in_nak / full_nak_deg
    balance_years = y * (1 - fraction)
    return lord_idx, balance_years

def generate_periods(start_date, lord_idx, total_years, level='dasa', max_depth=3):
    periods = []
    remaining = total_years
    i = lord_idx
    current_start = start_date
    total_cycle = 120
    depth_dict = {'dasa':0, 'bhukti':1, 'anthara':2, 'sukshma':3, 'prana':4, 'sub_prana':5}
    level_to_next = {0:'bhukti',1:'anthara',2:'sukshma',3:'prana',4:'sub_prana',5:None}
    depth = depth_dict.get(level, 0)
    while remaining > 0:
        lord = lords_full[i]
        lord_full = years[i]
        y = (lord_full / total_cycle) * total_years
        if y > remaining:
            y = remaining
        end_date = current_start + timedelta(days=y * 365.25)
        sub_periods = []
        if depth < max_depth - 1 and level_to_next.get(depth):
            next_level = level_to_next[depth]
            sub_periods = generate_periods(current_start, i, y, next_level, max_depth)
        periods.append((lord, current_start, end_date, sub_periods))
        remaining -= y
        current_start = end_date
        i = (i + 1) % 9
        if remaining <= 0:
            break
    return periods

def filter_from_birth(periods, birth_dt):
    filtered = []
    for lord, start, end, sub in periods:
        if end > birth_dt:
            adj_start = max(start, birth_dt)
            adj_sub = filter_from_birth(sub, birth_dt) if sub else []
            filtered.append((lord, adj_start, end, adj_sub))
    return filtered

def duration_str(delta, level='dasa'):
    total_days = delta.total_seconds() / 86400
    if total_days < 1 and level in ['sukshma', 'prana', 'sub_prana']:
        total_hours = total_days * 24
        hours = int(total_hours)
        minutes = int((total_hours - hours) * 60)
        if hours == 0 and minutes == 0:
            return "Less than 1 minute"
        return f"{hours}h {minutes}m"
    else:
        years = int(total_days / 365.25)
        rem_days = total_days % 365.25
        months = int(rem_days / 30.4375)
        days = int(rem_days % 30.4375)
        if years + months + days == 0:
            return "Less than 1 day"
        return f"{years}y {months}m {days}d"

def create_south_indian_chart(house_planets, lagna_sign, chart_type='Rasi'):
    """Create South Indian style chart using Streamlit columns"""
    st.markdown(f"<h4 style='text-align: center; color: #125336; margin-bottom: 20px;'>{chart_type} Chart - South Indian Style</h4>", unsafe_allow_html=True)
    
    # In South Indian style, houses are FIXED in position
    # Signs rotate based on lagna
    # Layout:
    #  12  1  2  3
    #  11  |  |  4
    #  10  |  |  5
    #   9  8  7  6
    
    sign_positions = {
        'Aries': 0, 'Taurus': 1, 'Gemini': 2, 'Cancer': 3,
        'Leo': 4, 'Virgo': 5, 'Libra': 6, 'Scorpio': 7,
        'Sagittarius': 8, 'Capricorn': 9, 'Aquarius': 10, 'Pisces': 11
    }
    
    lagna_idx = sign_positions[lagna_sign]
    
    # Create container with border
    st.markdown("""
    <style>
    .si-container {
        border: 3px solid #125336;
        padding: 0;
        background: white;
        max-width: 600px;
        margin: 0 auto;
    }
    .si-row {
        display: flex;
        width: 100%;
    }
    .si-cell {
        flex: 1;
        border: 1px solid #125336;
        padding: 10px;
        min-height: 120px;
        position: relative;
        background: white;
    }
    .si-center {
        background: #f8f9fa;
        border-left: 2px solid #125336;
        border-right: 2px solid #125336;
    }
    .si-house-num {
        position: absolute;
        top: 3px;
        left: 5px;
        font-size: 10px;
        color: #666;
        font-weight: bold;
    }
    .si-sign-name {
        position: absolute;
        top: 3px;
        right: 5px;
        font-size: 9px;
        color: #888;
    }
    .si-planets {
        margin-top: 25px;
        text-align: center;
        color: #125336;
        font-size: 11px;
        font-weight: 500;
    }
    </style>
    """, unsafe_allow_html=True)
    
    layout = [
        [12, 1, 2, 3],
        [11, 'X', 'X', 4],
        [10, 'X', 'X', 5],
        [9, 8, 7, 6]
    ]
    
    chart_html = '<div class="si-container">'
    
    for row_idx, row in enumerate(layout):
        chart_html += '<div class="si-row">'
        for col_idx, house_num in enumerate(row):
            if house_num == 'X':
                # Center cells with just borders
                chart_html += '<div class="si-cell si-center"></div>'
            else:
                # Calculate which sign is in this house
                sign_idx = (lagna_idx + house_num - 1) % 12
                sign = sign_names[sign_idx]
                
                # Get planets in this house
                planets = house_planets.get(house_num, [])
                planet_text = '<br>'.join(planets) if planets else 'Empty'
                
                chart_html += f'''
                <div class="si-cell">
                    <span class="si-house-num">{house_num}</span>
                    <span class="si-sign-name">{sign}</span>
                    <div class="si-planets">{planet_text}</div>
                </div>
                '''
        chart_html += '</div>'
    
    chart_html += '</div>'
    
    st.markdown(chart_html, unsafe_allow_html=True)

def create_north_indian_chart(house_planets, lagna_sign, chart_type='Rasi'):
    """Create North Indian style chart with diagonal lines"""
    st.markdown(f"<h4 style='text-align: center; color: #125336; margin-bottom: 20px;'>{chart_type} Chart - North Indian Style</h4>", unsafe_allow_html=True)
    
    # In North Indian chart, SIGNS are fixed in position
    # Houses rotate based on which sign is the lagna
    
    sign_positions = {
        'Aries': 0, 'Taurus': 1, 'Gemini': 2, 'Cancer': 3,
        'Leo': 4, 'Virgo': 5, 'Libra': 6, 'Scorpio': 7,
        'Sagittarius': 8, 'Capricorn': 9, 'Aquarius': 10, 'Pisces': 11
    }
    
    lagna_idx = sign_positions[lagna_sign]
    
    # North Indian diamond layout
    # Signs in fixed positions (clockwise from top):
    # Cancer, Leo, Virgo, Libra, Scorpio, Sagittarius, Capricorn, Aquarius, Pisces, Aries, Taurus, Gemini
    # Positions: Top, Top-Right, Right, Right-Bottom, Bottom, Bottom-Left, Left, Left-Top, and 4 corners
    
    st.markdown("""
    <style>
    .ni-container {
        width: 500px;
        height: 500px;
        margin: 0 auto;
        position: relative;
        background: white;
    }
    .ni-diamond {
        width: 100%;
        height: 100%;
        position: relative;
        transform: rotate(0deg);
    }
    /* Main square border */
    .ni-square {
        position: absolute;
        width: 350px;
        height: 350px;
        border: 3px solid #125336;
        top: 75px;
        left: 75px;
        transform: rotate(45deg);
        background: white;
    }
    /* Diagonal lines */
    .ni-line-v {
        position: absolute;
        width: 2px;
        height: 247px;
        background: #125336;
        top: 126px;
        left: 249px;
    }
    .ni-line-h {
        position: absolute;
        width: 247px;
        height: 2px;
        background: #125336;
        top: 249px;
        left: 126px;
    }
    /* House positions */
    .ni-house {
        position: absolute;
        padding: 5px;
        font-size: 10px;
        text-align: center;
        color: #125336;
    }
    .ni-house-num {
        font-weight: bold;
        font-size: 9px;
        color: #666;
    }
    .ni-sign {
        font-size: 8px;
        color: #888;
        margin-bottom: 3px;
    }
    .ni-planets {
        font-size: 9px;
        color: #125336;
        font-weight: 500;
        margin-top: 3px;
    }
    /* Specific house positions in diamond */
    .pos-top { top: 10px; left: 220px; width: 60px; }
    .pos-top-right { top: 50px; right: 80px; width: 70px; }
    .pos-right { top: 220px; right: 10px; width: 60px; }
    .pos-bottom-right { bottom: 80px; right: 50px; width: 70px; }
    .pos-bottom { bottom: 10px; left: 220px; width: 60px; }
    .pos-bottom-left { bottom: 80px; left: 50px; width: 70px; }
    .pos-left { top: 220px; left: 10px; width: 60px; }
    .pos-top-left { top: 50px; left: 80px; width: 70px; }
    .pos-top-inner { top: 100px; left: 195px; width: 110px; }
    .pos-right-inner { top: 195px; right: 100px; width: 110px; }
    .pos-bottom-inner { bottom: 100px; left: 195px; width: 110px; }
    .pos-left-inner { top: 195px; left: 100px; width: 110px; }
    .pos-center { top: 210px; left: 210px; width: 80px; background: #f0f7f4; border: 1px solid #125336; border-radius: 5px; font-weight: bold; }
    </style>
    """, unsafe_allow_html=True)
    
    # Fixed sign positions in North Indian chart (12 positions around diamond)
    sign_position_map = [
        ('Cancer', 'pos-top'),
        ('Leo', 'pos-top-right'), 
        ('Virgo', 'pos-right'),
        ('Libra', 'pos-bottom-right'),
        ('Scorpio', 'pos-bottom'),
        ('Sagittarius', 'pos-bottom-left'),
        ('Capricorn', 'pos-left'),
        ('Aquarius', 'pos-top-left'),
        ('Pisces', 'pos-top-inner'),
        ('Aries', 'pos-right-inner'),
        ('Taurus', 'pos-bottom-inner'),
        ('Gemini', 'pos-left-inner')
    ]
    
    chart_html = '<div class="ni-container"><div class="ni-diamond">'
    chart_html += '<div class="ni-square"></div>'
    chart_html += '<div class="ni-line-v"></div>'
    chart_html += '<div class="ni-line-h"></div>'
    
    # Place houses based on signs
    for sign_name, position_class in sign_position_map:
        sign_idx = sign_positions[sign_name]
        # Calculate which house this sign represents
        house_num = ((sign_idx - lagna_idx) % 12) + 1
        
        # Get planets in this house
        planets = house_planets.get(house_num, [])
        planet_text = '<br>'.join(planets[:3]) if planets else 'Empty'  # Limit to 3 planets for space
        
        chart_html += f'''
        <div class="ni-house {position_class}">
            <div class="ni-house-num">{house_num}</div>
            <div class="ni-sign">{sign_name[:3]}</div>
            <div class="ni-planets">{planet_text}</div>
        </div>
        '''
    
    # Add center ascendant marker
    chart_html += f'''
    <div class="ni-house pos-center">
        <div style="font-size: 11px;">Asc</div>
        <div style="font-size: 10px;">{lagna_sign}</div>
    </div>
    '''
    
    chart_html += '</div></div>'
    
    st.markdown(chart_html, unsafe_allow_html=True)

def compute_chart(name, date_obj, time_str, lat, lon, tz_offset, max_depth):
    # Parse time string HH:MM
    try:
        time_parts = time_str.split(':')
        hour = int(time_parts[0])
        minute = int(time_parts[1]) if len(time_parts) > 1 else 0
        if hour < 0 or hour > 23 or minute < 0 or minute > 59:
            raise ValueError("Invalid time format")
    except:
        raise ValueError("Time must be in HH:MM format (24-hour)")
    
    local_dt = datetime.combine(date_obj, datetime.min.time().replace(hour=hour, minute=minute))
    utc_dt = local_dt - timedelta(hours=tz_offset)
    t = Time(utc_dt)
    jd = t.jd
    year = utc_dt.year
    ayan = get_lahiri_ayanamsa(year)

    with solar_system_ephemeris.set('builtin'):
        planet_bodies = {
            'sun': 'sun',
            'moon': 'moon',
            'mercury': 'mercury',
            'venus': 'venus',
            'mars': 'mars',
            'jupiter': 'jupiter',
            'saturn': 'saturn'
        }
        lon_trop = {}
        for name_p, body in planet_bodies.items():
            p = get_body(body, t)
            ecl = p.transform_to(GeocentricTrueEcliptic())
            lon_trop[name_p] = ecl.lon.deg

    d = jd - 2451545.0
    T = d / 36525.0
    omega = (125.04452 - 1934.136261 * T + 0.0020708 * T**2 + T**3 / 450000) % 360
    lon_trop['rahu'] = omega
    lon_trop['ketu'] = (omega + 180) % 360

    lon_sid = {p: get_sidereal_lon(lon_trop[p], ayan) for p in lon_trop}

    asc_trop = get_ascendant(jd, lat, lon)
    lagna_sid = get_sidereal_lon(asc_trop, ayan)

    data = []
    asc_deg = lagna_sid % 360
    asc_sign = get_sign(asc_deg)
    asc_nak, asc_pada, asc_ld, asc_sl = get_nakshatra_details(asc_deg)
    data.append(['Asc', f"{asc_deg:.2f}", asc_sign, asc_nak, asc_pada, f"{asc_ld}/{asc_sl}"])

    for p in ['sun', 'moon', 'mars', 'mercury', 'jupiter', 'venus', 'saturn', 'rahu', 'ketu']:
        p_name = p.capitalize()
        lon = lon_sid[p]
        sign = get_sign(lon)
        nak, pada, ld, sl = get_nakshatra_details(lon)
        data.append([p_name, f"{lon:.2f}", sign, nak, pada, f"{ld}/{sl}"])

    df_planets = pd.DataFrame(data, columns=['Planet', 'Deg', 'Sign', 'Nakshatra', 'Pada', 'Ld/SL'])

    house_planets_rasi = defaultdict(list)
    positions = {**lon_sid, 'asc': lagna_sid}
    for p, lon in positions.items():
        h = get_house(lon, lagna_sid)
        house_planets_rasi[h].append(p.capitalize() if p != 'asc' else 'Asc')

    rasi_data = []
    for h in range(1, 13):
        sign_start = (lagna_sid + (h - 1) * 30) % 360
        sign = get_sign(sign_start)
        pls = ', '.join(sorted(house_planets_rasi[h])) if house_planets_rasi[h] else 'Empty'
        rasi_data.append([f"House {h}", sign, pls])
    df_rasi = pd.DataFrame(rasi_data, columns=['House', 'Sign', 'Planets'])

    nav_lagna = (lagna_sid * 9) % 360
    house_planets_nav = defaultdict(list)
    for p, lon in lon_sid.items():
        nav_lon = (lon * 9) % 360
        nav_s = int(nav_lon / 30)
        nav_l = int(nav_lagna / 30)
        nav_h = (nav_s - nav_l) % 12 + 1
        house_planets_nav[nav_h].append(p.capitalize())

    nav_data = []
    for h in range(1, 13):
        nav_sign_start = (nav_lagna + (h - 1) * 30) % 360
        nav_sign = get_sign(nav_sign_start)
        pls = ', '.join(sorted(house_planets_nav[h])) if house_planets_nav[h] else 'Empty'
        nav_data.append([f"House {h}", nav_sign, pls])
    df_nav = pd.DataFrame(nav_data, columns=['House', 'Sign', 'Planets'])

    lagna_sign = get_sign(lagna_sid)
    lagna_idx = sign_names.index(lagna_sign)
    planet_to_house = {p.capitalize(): get_house(lon_sid[p], lagna_sid) for p in lon_sid}
    aspects_dict = {
        'Sun': [7],
        'Moon': [7],
        'Mars': [4,7,8],
        'Mercury': [7],
        'Jupiter': [5,7,9],
        'Venus': [7],
        'Saturn': [3,7,10],
    }
    house_status_data = []
    for h in range(1,13):
        sign_idx = (lagna_idx + h - 1) % 12
        lord = sign_lords[sign_idx]
        lord_house = planet_to_house[lord]
        aspecting = []
        for planet, offsets in aspects_dict.items():
            if planet in planet_to_house:
                p_h = planet_to_house[planet]
                for off in offsets:
                    a_h = ((p_h - 1 + (off - 1)) % 12) + 1
                    if a_h == h:
                        aspecting.append(planet)
        aspect_str = ', '.join(aspecting) if aspecting else 'None'
        pls = ', '.join(sorted(house_planets_rasi[h])) if house_planets_rasi[h] else 'Empty'
        house_status_data.append([f"House {h}", pls, aspect_str, lord, f"House {lord_house}"])
    df_house_status = pd.DataFrame(house_status_data, columns=['House', 'Planets', 'Aspects from', 'Lord', 'Lord in'])

    moon_lon = lon_sid['moon']
    dasa_start_idx, balance_years = generate_vimshottari_dasa(moon_lon)
    full_first = years[dasa_start_idx]
    passed_years = full_first - balance_years
    dasa_start_dt = utc_dt - timedelta(days=passed_years * 365.25)

    dasa_periods = generate_periods(dasa_start_dt, dasa_start_idx, 120, level='dasa', max_depth=max_depth)
    dasa_periods_filtered = filter_from_birth(dasa_periods, utc_dt)

    max_depth_options = {1: 'Dasa only', 2: 'Dasa + Bhukti', 3: 'Dasa + Bhukti + Anthara', 
                         4: 'Dasa + Bhukti + Anthara + Sukshma', 5: 'Dasa + Bhukti + Anthara + Sukshma + Prana',
                         6: 'Dasa + Bhukti + Anthara + Sukshma + Prana + Sub-Prana'}
    selected_depth = max_depth_options[max_depth]

    nav_lagna_sign = get_sign(nav_lagna)
    
    # Moon nakshatra for summary
    moon_nak, moon_pada, moon_ld, moon_sl = get_nakshatra_details(moon_lon)

    return {
        'name': name,
        'df_planets': df_planets,
        'df_rasi': df_rasi,
        'df_nav': df_nav,
        'df_house_status': df_house_status,
        'dasa_periods_filtered': dasa_periods_filtered,
        'lagna_sid': lagna_sid,
        'nav_lagna': nav_lagna,
        'lagna_sign': lagna_sign,
        'nav_lagna_sign': nav_lagna_sign,
        'moon_rasi': get_sign(moon_lon),
        'moon_nakshatra': moon_nak,
        'moon_pada': moon_pada,
        'house_planets_rasi': dict(house_planets_rasi),
        'house_planets_nav': dict(house_planets_nav),
        'selected_depth': selected_depth,
        'utc_dt': utc_dt,
        'max_depth': max_depth
    }

# Streamlit UI
st.set_page_config(page_title="Vedic Astrology", layout="wide")

# Enhanced CSS
st.markdown("""
<style>
    .stApp {
        background-color: white;
        color: #125336;
    }
    .stTextInput > div > div > input,
    .stSelectbox > div > div > select,
    .stNumberInput > div > div > input {
        background-color: white;
        color: #125336;
        border: 1px solid #125336;
    }
    .stButton > button {
        background-color: #125336;
        color: white;
        border: none;
        padding: 0.5rem 2rem;
        font-size: 1.1rem;
        font-weight: 600;
    }
    .stButton > button:hover {
        background-color: #0a3d22;
        box-shadow: 0 4px 6px rgba(0,0,0,0.1);
    }
    h1, h2, h3 {
        color: #125336 !important;
    }
    .summary-box {
        background-color: #f0f7f4;
        padding: 1.5rem;
        border-radius: 10px;
        border: 2px solid #125336;
        margin: 1.5rem 0;
    }
    .summary-box h3 {
        margin-top: 0;
        color: #125336 !important;
    }
    .summary-item {
        font-size: 1.1rem;
        margin: 0.5rem 0;
        color: #125336;
    }
    div[data-testid="stVerticalBlock"] > div:first-child {
        padding-top: 0 !important;
    }
</style>
""", unsafe_allow_html=True)

st.title("Vedic Astrology Chart Generator")

# Initialize session state
if 'chart_data' not in st.session_state:
    st.session_state.chart_data = None
if 'location_data' not in st.session_state:
    st.session_state.location_data = None
if 'search_results' not in st.session_state:
    st.session_state.search_results = []
if 'chart_style' not in st.session_state:
    st.session_state.chart_style = 'South Indian'

# Initialize geolocator
@st.cache_resource
def get_geolocator():
    geolocator = Nominatim(user_agent="vedic_astro_app")
    return RateLimiter(geolocator.geocode, min_delay_seconds=1)

geocode = get_geolocator()

# Input Section
st.subheader("Birth Details")

# Name input
name = st.text_input("Name", placeholder="Enter full name")

col1, col2, col3 = st.columns(3)

with col1:
    birth_date = st.date_input(
        "Birth Date",
        value=datetime.now().date(),
        min_value=datetime(1900, 1, 1).date(),
        max_value=datetime.now().date()
    )

with col2:
    birth_time = st.text_input(
        "Birth Time (HH:MM in 24-hour format)",
        placeholder="14:30",
        help="Example: 14:30 for 2:30 PM"
    )

with col3:
    tz_offset = st.number_input(
        "Timezone offset (hrs)",
        value=5.5,
        step=0.5,
        help="Offset from UTC (e.g., IST = 5.5)"
    )

# Custom coordinates checkbox
use_custom_coords = st.checkbox("Custom latitude and longitude?")

if use_custom_coords:
    col_lat, col_lon = st.columns(2)
    with col_lat:
        lat = st.number_input("Latitude", value=13.08, format="%.4f")
    with col_lon:
        lon = st.number_input("Longitude", value=80.27, format="%.4f")
    location_data = {'lat': lat, 'lon': lon}
else:
    # City search with autocomplete
    city_query = st.text_input(
        "Search City",
        placeholder="Start typing city name...",
        key="city_input"
    )
    
    # Automatically search as user types
    if city_query and len(city_query) >= 2:
        try:
            locations = geocode(city_query, exactly_one=False, limit=5)
            if locations:
                st.session_state.search_results = [
                    {
                        'display': f"{loc.address}",
                        'lat': loc.latitude,
                        'lon': loc.longitude,
                        'address': loc.address
                    }
                    for loc in locations
                ]
            else:
                city_key = city_query.title()
                if city_key in cities_fallback:
                    st.session_state.search_results = [{
                        'display': f"{city_key} (Fallback)",
                        'lat': cities_fallback[city_key]['lat'],
                        'lon': cities_fallback[city_key]['lon'],
                        'address': city_key
                    }]
                else:
                    st.session_state.search_results = []
        except:
            city_key = city_query.title()
            if city_key in cities_fallback:
                st.session_state.search_results = [{
                    'display': f"{city_key} (Fallback)",
                    'lat': cities_fallback[city_key]['lat'],
                    'lon': cities_fallback[city_key]['lon'],
                    'address': city_key
                }]
            else:
                st.session_state.search_results = []
    elif len(city_query) < 2:
        st.session_state.search_results = []
    
    # Show search results
    if st.session_state.search_results:
        result_options = [r['display'] for r in st.session_state.search_results]
        selected = st.selectbox(
            "Select location",
            options=result_options,
            key="location_selector"
        )
        selected_idx = result_options.index(selected)
        location_data = st.session_state.search_results[selected_idx]
        lat = location_data['lat']
        lon = location_data['lon']
        st.success(f"Selected: {location_data['address']} (Lat: {lat:.2f}, Lon: {lon:.2f})")
    else:
        lat, lon = 13.08, 80.27
        location_data = {'lat': lat, 'lon': lon}
        st.info("Using default location: Chennai, India")

# Dasa depth selector
col_depth, col_style = st.columns(2)
with col_depth:
    max_depth_options = {
        1: 'Dasa only',
        2: 'Dasa + Bhukti',
        3: 'Dasa + Bhukti + Anthara',
        4: 'Dasa + Bhukti + Anthara + Sukshma',
        5: 'Dasa + Bhukti + Anthara + Sukshma + Prana',
        6: 'Dasa + Bhukti + Anthara + Sukshma + Prana + Sub-Prana'
    }
    selected_depth_str = st.selectbox(
        "Period Depth",
        options=list(max_depth_options.values()),
        index=2
    )
    max_depth = list(max_depth_options.keys())[list(max_depth_options.values()).index(selected_depth_str)]

with col_style:
    chart_style = st.selectbox(
        "Chart Style",
        options=["South Indian", "North Indian"],
        index=0
    )

# Generate button
if st.button("Generate Chart", use_container_width=True):
    if not name:
        st.error("Please enter a name.")
    elif not birth_time:
        st.error("Please enter birth time in HH:MM format.")
    else:
        try:
            with st.spinner("Calculating chart..."):
                st.session_state.chart_data = compute_chart(
                    name, birth_date, birth_time, lat, lon, tz_offset, max_depth
                )
                st.session_state.chart_style = chart_style
            st.success("Chart generated successfully!")
            st.rerun()
        except ValueError as e:
            st.error(f"Invalid input: {e}")
        except Exception as e:
            st.error(f"Error generating chart: {e}")

# Display results
if st.session_state.chart_data:
    chart_data = st.session_state.chart_data
    chart_style = st.session_state.get('chart_style', 'South Indian')
    
    st.markdown("---")
    
    # Summary Box - First thing displayed
    st.markdown(f"""
    <div class="summary-box">
        <h3>Chart Summary</h3>
        <div class="summary-item"><strong>Name:</strong> {chart_data['name']}</div>
        <div class="summary-item"><strong>Lagna:</strong> {chart_data['lagna_sign']} ({chart_data['lagna_sid']:.2f}°)</div>
        <div class="summary-item"><strong>Rasi (Moon Sign):</strong> {chart_data['moon_rasi']}</div>
        <div class="summary-item"><strong>Nakshatra:</strong> {chart_data['moon_nakshatra']} (Pada {chart_data['moon_pada']})</div>
    </div>
    """, unsafe_allow_html=True)
    
    # Planetary Details
    st.subheader("Planetary Positions")
    st.dataframe(chart_data['df_planets'], hide_index=True, use_container_width=True)
    
    # Rasi Chart - Visual
    st.subheader("Rasi Chart (D1)")
    if chart_style == "South Indian":
        create_south_indian_chart(chart_data['house_planets_rasi'], chart_data['lagna_sign'], 'Rasi')
    else:
        create_north_indian_chart(chart_data['house_planets_rasi'], chart_data['lagna_sign'], 'Rasi')
    
    with st.expander("View Rasi Chart as Table"):
        st.dataframe(chart_data['df_rasi'], hide_index=True, use_container_width=True)
    
    # House Status
    st.subheader("House Analysis")
    st.dataframe(chart_data['df_house_status'], hide_index=True, use_container_width=True)
    
    # Navamsa Chart - Visual
    st.subheader("Navamsa Chart (D9)")
    if chart_style == "South Indian":
        create_south_indian_chart(chart_data['house_planets_nav'], chart_data['nav_lagna_sign'], 'Navamsa')
    else:
        create_north_indian_chart(chart_data['house_planets_nav'], chart_data['nav_lagna_sign'], 'Navamsa')
    
    with st.expander("View Navamsa Chart as Table"):
        st.write(f"Navamsa Lagna: {chart_data['nav_lagna_sign']} ({chart_data['nav_lagna']:.2f}°)")
        st.dataframe(chart_data['df_nav'], hide_index=True, use_container_width=True)
    
    # Vimshottari Dasa
    st.subheader(f"Vimshottari Dasa ({chart_data['selected_depth']})")
    
    dasa_periods_filtered = chart_data['dasa_periods_filtered']
    utc_dt = chart_data['utc_dt']
    max_depth = chart_data['max_depth']
    
    # Main Dasa table
    dasa_data = []
    for lord, start, end, _ in dasa_periods_filtered:
        dur = end - start
        dasa_data.append({
            'Planet': lord,
            'Start': start.strftime('%Y-%m-%d'),
            'End': end.strftime('%Y-%m-%d'),
            'Duration': duration_str(dur, 'dasa')
        })
    df_dasa = pd.DataFrame(dasa_data)
    st.dataframe(df_dasa, hide_index=True, use_container_width=True)
    
    # Nested periods with expanders
    if max_depth >= 2:
        with st.expander("View Bhuktis (Sub-periods)", expanded=False):
            if dasa_periods_filtered:
                dasa_options = [f"{p[0]} ({p[1].strftime('%Y-%m-%d')} - {p[2].strftime('%Y-%m-%d')})" 
                               for p in dasa_periods_filtered]
                selected_dasa = st.selectbox("Select Dasa:", dasa_options, key="dasa_select")
                sel_idx = dasa_options.index(selected_dasa)
                bhuktis = dasa_periods_filtered[sel_idx][3]
                
                bhukti_data = []
                for b_lord, b_start, b_end, _ in bhuktis:
                    dur = b_end - b_start
                    bhukti_data.append({
                        'Planet': b_lord,
                        'Start': b_start.strftime('%Y-%m-%d'),
                        'End': b_end.strftime('%Y-%m-%d'),
                        'Duration': duration_str(dur, 'bhukti')
                    })
                df_bhukti = pd.DataFrame(bhukti_data)
                st.dataframe(df_bhukti, hide_index=True, use_container_width=True)
                
                if max_depth >= 3:
                    with st.expander("View Antharas", expanded=False):
                        if bhuktis:
                            bhukti_options = [f"{p[0]} ({p[1].strftime('%Y-%m-%d')} - {p[2].strftime('%Y-%m-%d')})" 
                                            for p in bhuktis]
                            selected_bhukti = st.selectbox("Select Bhukti:", bhukti_options, key="bhukti_select")
                            sel_b_idx = bhukti_options.index(selected_bhukti)
                            antharas = bhuktis[sel_b_idx][3]
                            
                            anthara_data = []
                            for a_lord, a_start, a_end, _ in antharas:
                                dur = a_end - a_start
                                anthara_data.append({
                                    'Planet': a_lord,
                                    'Start': a_start.strftime('%Y-%m-%d %H:%M'),
                                    'End': a_end.strftime('%Y-%m-%d %H:%M'),
                                    'Duration': duration_str(dur, 'anthara')
                                })
                            df_anthara = pd.DataFrame(anthara_data)
                            st.dataframe(df_anthara, hide_index=True, use_container_width=True)
                            
                            if max_depth >= 4:
                                with st.expander("View Sukshmas", expanded=False):
                                    if antharas:
                                        anthara_options = [f"{p[0]} ({p[1].strftime('%Y-%m-%d %H:%M')} - {p[2].strftime('%Y-%m-%d %H:%M')})" 
                                                         for p in antharas]
                                        selected_anthara = st.selectbox("Select Anthara:", anthara_options, key="anthara_select")
                                        sel_a_idx = anthara_options.index(selected_anthara)
                                        sukshmas = antharas[sel_a_idx][3]
                                        
                                        sukshma_data = []
                                        for s_lord, s_start, s_end, _ in sukshmas:
                                            dur = s_end - s_start
                                            sukshma_data.append({
                                                'Planet': s_lord,
                                                'Start': s_start.strftime('%Y-%m-%d %H:%M'),
                                                'End': s_end.strftime('%Y-%m-%d %H:%M'),
                                                'Duration': duration_str(dur, 'sukshma')
                                            })
                                        df_sukshma = pd.DataFrame(sukshma_data)
                                        st.dataframe(df_sukshma, hide_index=True, use_container_width=True)
                                        
                                        if max_depth >= 5:
                                            with st.expander("View Pranas", expanded=False):
                                                if sukshmas:
                                                    sukshma_options = [f"{p[0]} ({p[1].strftime('%Y-%m-%d %H:%M')} - {p[2].strftime('%Y-%m-%d %H:%M')})" 
                                                                     for p in sukshmas]
                                                    selected_sukshma = st.selectbox("Select Sukshma:", sukshma_options, key="sukshma_select")
                                                    sel_s_idx = sukshma_options.index(selected_sukshma)
                                                    pranas = sukshmas[sel_s_idx][3]
                                                    
                                                    prana_data = []
                                                    for pr_lord, pr_start, pr_end, _ in pranas:
                                                        dur = pr_end - pr_start
                                                        prana_data.append({
                                                            'Planet': pr_lord,
                                                            'Start': pr_start.strftime('%Y-%m-%d %H:%M'),
                                                            'End': pr_end.strftime('%Y-%m-%d %H:%M'),
                                                            'Duration': duration_str(dur, 'prana')
                                                        })
                                                    df_prana = pd.DataFrame(prana_data)
                                                    st.dataframe(df_prana, hide_index=True, use_container_width=True)
                                                    
                                                    if max_depth >= 6:
                                                        with st.expander("View Sub-Pranas", expanded=False):
                                                            if pranas:
                                                                prana_options = [f"{p[0]} ({p[1].strftime('%Y-%m-%d %H:%M')} - {p[2].strftime('%Y-%m-%d %H:%M')})" 
                                                                               for p in pranas]
                                                                selected_prana = st.selectbox("Select Prana:", prana_options, key="prana_select")
                                                                sel_pr_idx = prana_options.index(selected_prana)
                                                                sub_pranas = pranas[sel_pr_idx][3]
                                                                
                                                                sub_prana_data = []
                                                                for sp_lord, sp_start, sp_end, _ in sub_pranas:
                                                                    dur = sp_end - sp_start
                                                                    sub_prana_data.append({
                                                                        'Planet': sp_lord,
                                                                        'Start': sp_start.strftime('%Y-%m-%d %H:%M'),
                                                                        'End': sp_end.strftime('%Y-%m-%d %H:%M'),
                                                                        'Duration': duration_str(dur, 'sub_prana')
                                                                    })
                                                                df_sub_prana = pd.DataFrame(sub_prana_data)
                                                                st.dataframe(df_sub_prana, hide_index=True, use_container_width=True)
    
    st.info("Note: Periods filtered from birth. Durations approximate.")
else:
    st.info("Enter details above and click 'Generate Chart' to begin")

st.markdown("---")
st.caption("Vedic Astrology Chart Generator")
