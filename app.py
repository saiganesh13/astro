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
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch
import io  # for optional SVG rendering

# Note: pip install streamlit astropy geopy pandas matplotlib

# -------- MATPLOTLIB DEFAULTS (CRISP & THIN) --------
plt.rcParams.update({
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "lines.linewidth": 0.28
})

# -------- CONSTANTS --------
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
nak_names = [
    'Ashwini','Bharani','Krittika','Rohini','Mrigashira','Ardra','Punarvasu','Pushya','Ashlesha',
    'Magha','Purva Phalguni','Uttara Phalguni','Hasta','Chitra','Swati','Vishakha','Anuradha',
    'Jyeshta','Mula','Purva Ashadha','Uttara Ashadha','Shravana','Dhanishta','Shatabhisha',
    'Purva Bhadrapada','Uttara Bhadrapada','Revati'
]
years = [7, 20, 6, 10, 7, 18, 16, 19, 17] * 3
sign_lords = ['Mars','Venus','Mercury','Moon','Sun','Mercury','Venus','Mars','Jupiter','Saturn','Saturn','Jupiter']

# -------- ASTRONOMY/AYANAMSA UTILS --------
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
        bodies = {'sun':'sun','moon':'moon','mercury':'mercury','venus':'venus','mars':'mars','jupiter':'jupiter','saturn':'saturn'}
        lon_trop = {}
        for nm, body in bodies.items():
            p = get_body(body, t)
            ecl = p.transform_to(GeocentricTrueEcliptic())
            lon_trop[nm] = ecl.lon.deg

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

    for p in ['sun','moon','mars','mercury','jupiter','venus','saturn','rahu','ketu']:
        p_name = p.capitalize()
        lon = lon_sid[p]
        sign = get_sign(lon)
        nak, pada, ld, sl = get_nakshatra_details(lon)
        data.append([p_name, f"{lon:.2f}", sign, nak, pada, f"{ld}/{sl}"])

    df_planets = pd.DataFrame(data, columns=['Planet','Deg','Sign','Nakshatra','Pada','Ld/SL'])

    house_planets_rasi = defaultdict(list)
    positions = {**lon_sid, 'asc': lagna_sid}
    for p, lon in positions.items():
        h = get_house(lon, lagna_sid)
        house_planets_rasi[h].append(p.capitalize() if p != 'asc' else 'Asc')

    rasi_data = []
    for h in range(1, 13):
        sign_start = (lagna_sid + (h - 1) * 30) % 360
        rasi_data.append([f"House {h}", get_sign(sign_start),
                          ', '.join(sorted(house_planets_rasi[h])) if house_planets_rasi[h] else 'Empty'])
    df_rasi = pd.DataFrame(rasi_data, columns=['House','Sign','Planets'])

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
        nav_data.append([f"House {h}", get_sign(nav_sign_start),
                         ', '.join(sorted(house_planets_nav[h])) if house_planets_nav[h] else 'Empty'])
    df_nav = pd.DataFrame(nav_data, columns=['House','Sign','Planets'])

    lagna_sign = get_sign(lagna_sid)
    lagna_idx = sign_names.index(lagna_sign)
    planet_to_house = {p.capitalize(): get_house(lon_sid[p], lagna_sid) for p in lon_sid}
    aspects_dict = {'Sun':[7],'Moon':[7],'Mars':[4,7,8],'Mercury':[7],'Jupiter':[5,7,9],'Venus':[7],'Saturn':[3,7,10]}
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
    df_house_status = pd.DataFrame(house_status_data, columns=['House','Planets','Aspects from','Lord','Lord in'])

    moon_lon = lon_sid['moon']
    dasa_start_idx, balance_years = generate_vimshottari_dasa(moon_lon)
    full_first = years[dasa_start_idx]
    passed_years = full_first - balance_years
    dasa_start_dt = utc_dt - timedelta(days=passed_years * 365.25)
    dasa_periods = generate_periods(dasa_start_dt, dasa_start_idx, 120, level='dasa', max_depth=max_depth)
    dasa_periods_filtered = filter_from_birth(dasa_periods, utc_dt)

    max_depth_options = {
        1:'Dasa only', 2:'Dasa + Bhukti', 3:'Dasa + Bhukti + Anthara',
        4:'Dasa + Bhukti + Anthara + Sukshma', 5:'Dasa + Bhukti + Anthara + Sukshma + Prana',
        6:'Dasa + Bhukti + Anthara + Sukshma + Prana + Sub-Prana'
    }

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
        'nav_lagna_sign': get_sign(nav_lagna),
        'moon_rasi': get_sign(moon_lon),
        'moon_nakshatra': get_nakshatra_details(moon_lon)[0],
        'moon_pada': get_nakshatra_details(moon_lon)[1],
        'selected_depth': max_depth_options[max_depth],
        'utc_dt': utc_dt,
        'max_depth': max_depth,
        'house_to_planets_rasi': house_planets_rasi,
        'house_to_planets_nav': house_planets_nav
    }

# -------- PLOTTERS (plain titles, non-bold planets, left-top sign, dynamic spacing) --------
def plot_north_indian_style(ax, house_to_planets, house_to_sign, title):
    scale = 0.8
    house_positions = {
        1:(0,0.5*scale),2:(-0.5*scale,0.25*scale),3:(-0.75*scale,0),4:(-0.5*scale,-0.25*scale),
        5:(0,-0.5*scale),6:(0.5*scale,-0.25*scale),7:(0.75*scale,0),8:(0.5*scale,0.25*scale),
        9:(0.25*scale,0.5*scale),10:(0.5*scale,0.75*scale),11:(-0.25*scale,0.75*scale),12:(-0.5*scale,0.5*scale)
    }
    ax.add_patch(patches.RegularPolygon((0,0),4,radius=0.8*scale,orientation=radians(45),
                                        edgecolor='black',facecolor='none',linewidth=0.3))
    for X1,Y1,X2,Y2 in [(0,-0.8*scale,0,0.8*scale),(-0.8*scale,0,0.8*scale,0),
                        (-0.4*scale,0.4*scale,0.4*scale,0.4*scale),(-0.4*scale,-0.4*scale,0.4*scale,-0.4*scale)]:
        ax.plot([X1,X2],[Y1,Y2],'k-',linewidth=0.22)
    ax.plot([-0.4*scale,-0.4*scale],[-0.4*scale,0.4*scale],'k-',linewidth=0.22)
    ax.plot([0.4*scale,0.4*scale],[-0.4*scale,0.4*scale],'k-',linewidth=0.22)

    box_size, half = 0.16, 0.08
    pad = 0.01

    for h in range(1,13):
        x,y = house_positions[h]
        sign = house_to_sign.get(h,'')
        planets_list = sorted(house_to_planets.get(h,[]))

        ax.add_patch(FancyBboxPatch((x-half,y-half),box_size,box_size,
                                    boxstyle="round,pad=0.004",ec="black",fc="#F5F5F5",
                                    alpha=0.88,linewidth=0.3))

        # left-top sign
        ax.text(x-half+pad, y+half-pad, sign[:3], ha='left', va='top', fontsize=3.2)

        if planets_list:
            avail = box_size - (pad + 0.028)       # space under sign
            n = len(planets_list)
            line_h = min(0.019, max(0.012, avail / max(n,1)))
            start_y = y - half + pad + 0.02
            for i, planet in enumerate(planets_list):
                py = start_y + i*line_h
                ax.text(x, py, planet[:5], ha='center', va='center', fontsize=3.6)  # plain text

    ax.set_xlim(-1,1); ax.set_ylim(-1,1); ax.set_aspect('equal')
    ax.set_title(title, fontsize=4.6, fontweight='normal')
    ax.axis('off')

def plot_south_indian_style(ax, house_to_planets, lagna_sign, title):
    sign_positions = {'Pisces':(0,3),'Aries':(1,3),'Taurus':(2,3),'Gemini':(3,3),
                      'Cancer':(3,2),'Leo':(3,1),'Virgo':(3,0),
                      'Libra':(2,0),'Scorpio':(1,0),'Sagittarius':(0,0),
                      'Capricorn':(0,1),'Aquarius':(0,2)}
    lagna_idx = sign_names.index(lagna_sign)
    house_for_sign = {s: ((i - lagna_idx) % 12) + 1 for i, s in enumerate(sign_names)}

    box_w, box_h, spacing = 0.46, 0.46, 0.52
    pad = 0.02
    for sign,(gx,gy) in sign_positions.items():
        h = house_for_sign[sign]
        planets_list = sorted(house_to_planets.get(h,[]))
        x = gx*spacing + 0.22
        y = (3-gy)*spacing + 0.22

        ax.add_patch(FancyBboxPatch((x,y),box_w,box_h,boxstyle="round,pad=0.004",
                                    ec="black",fc="#F5F5F5",alpha=0.92,linewidth=0.32))
        ax.text(x+pad, y+pad, sign[:3], ha='left', va='top', fontsize=3.2)

        if planets_list:
            avail = box_h - (pad + 0.04)
            n = len(planets_list)
            line_h = min(0.045, max(0.028, avail / max(n,1)))
            start_y = y + pad + 0.035
            for i, planet in enumerate(planets_list):
                py = start_y + i*line_h
                ax.text(x + box_w/2, py, planet, ha='center', va='top', fontsize=3.6)  # plain text

    ax.set_xlim(0,3); ax.set_ylim(0,3); ax.set_aspect('equal'); ax.invert_yaxis()
    ax.set_title(title, fontsize=4.6, fontweight='normal')
    ax.axis('off')

# -------- STREAMLIT UI --------
st.set_page_config(page_title="Sivapathy Horoscope", layout="wide")

st.markdown("""
<style>
    .stApp { background-color: white; color: #125336; }
    .stTextInput > div > div > input,
    .stSelectbox > div > div > select,
    .stNumberInput > div > div > input {
        background-color: white; color: #125336; border: 1px solid #125336;
    }
    .stButton > button {
        background-color: #125336; color: white; border: none; padding: 0.5rem 2rem;
        font-size: 1.1rem; font-weight: 600;
    }
    .stButton > button:hover { background-color: #0a3d22; box-shadow: 0 4px 6px rgba(0,0,0,0.1); }
    h1, h2, h3 { color: #125336 !important; }
    .summary-box { background-color: #f0f7f4; padding: 1.2rem; border-radius: 10px; border: 2px solid #125336; margin: 1rem 0; }
    .summary-item { font-size: 1.05rem; margin: 0.35rem 0; color: #125336; }
    div[data-testid="stVerticalBlock"] > div:first-child { padding-top: 0 !important; }
</style>
""", unsafe_allow_html=True)

st.title("Sivapathy Astrology Data Generator")

# Session
if 'chart_data' not in st.session_state: st.session_state.chart_data = None
if 'location_data' not in st.session_state: st.session_state.location_data = None
if 'search_results' not in st.session_state: st.session_state.search_results = []

@st.cache_resource
def get_geolocator():
    geolocator = Nominatim(user_agent="vedic_astro_app")
    return RateLimiter(geolocator.geocode, min_delay_seconds=1)

geocode = get_geolocator()

# Inputs
st.subheader("Birth Details")
name = st.text_input("Name", placeholder="Enter full name")
c1, c2, c3 = st.columns(3)
with c1:
    birth_date = st.date_input("Birth Date", value=datetime.now().date(),
                               min_value=datetime(1900,1,1).date(), max_value=datetime.now().date())
with c2:
    birth_time = st.text_input("Birth Time (HH:MM in 24-hour format)",
                               placeholder="14:30", help="Example: 14:30 for 2:30 PM")
with c3:
    tz_offset = st.number_input("Timezone offset (hrs)", value=5.5, step=0.5,
                                help="Offset from UTC (e.g., IST = 5.5)")

use_custom_coords = st.checkbox("Custom latitude and longitude?")
if use_custom_coords:
    clat, clon = st.columns(2)
    with clat: lat = st.number_input("Latitude", value=13.08, format="%.4f")
    with clon: lon = st.number_input("Longitude", value=80.27, format="%.4f")
    location_data = {'lat': lat, 'lon': lon}
else:
    city_query = st.text_input("Search City", placeholder="Start typing city name...", key="city_input")
    if city_query and len(city_query) >= 2:
        try:
            locations = geocode(city_query, exactly_one=False, limit=5)
            if locations:
                st.session_state.search_results = [{'display': f"{loc.address}", 'lat': loc.latitude, 'lon': loc.longitude, 'address': loc.address} for loc in locations]
            else:
                city_key = city_query.title()
                if city_key in cities_fallback:
                    st.session_state.search_results = [{'display': f"{city_key} (Fallback)",
                                                        'lat': cities_fallback[city_key]['lat'],
                                                        'lon': cities_fallback[city_key]['lon'],
                                                        'address': city_key}]
                else:
                    st.session_state.search_results = []
        except:
            city_key = city_query.title()
            if city_key in cities_fallback:
                st.session_state.search_results = [{'display': f"{city_key} (Fallback)",
                                                    'lat': cities_fallback[city_key]['lat'],
                                                    'lon': cities_fallback[city_key]['lon'],
                                                    'address': city_key}]
            else:
                st.session_state.search_results = []
    elif len(city_query) < 2:
        st.session_state.search_results = []

    if st.session_state.search_results:
        options = [r['display'] for r in st.session_state.search_results]
        selected = st.selectbox("Select location", options=options, key="location_selector")
        idx = options.index(selected)
        location_data = st.session_state.search_results[idx]
        lat = location_data['lat']; lon = location_data['lon']
        st.success(f"Selected: {location_data['address']} (Lat: {lat:.2f}, Lon: {lon:.2f})")
    else:
        lat, lon = 13.08, 80.27
        location_data = {'lat': lat, 'lon': lon}
        st.info("Using default location: Chennai, India")

max_depth_options = {
    1:'Dasa only', 2:'Dasa + Bhukti', 3:'Dasa + Bhukti + Anthara',
    4:'Dasa + Bhukti + Anthara + Sukshma', 5:'Dasa + Bhukti + Anthara + Sukshma + Prana',
    6:'Dasa + Bhukti + Anthara + Sukshma + Prana + Sub-Prana'
}
selected_depth_str = st.selectbox("Period Depth", options=list(max_depth_options.values()), index=2)
max_depth = list(max_depth_options.keys())[list(max_depth_options.values()).index(selected_depth_str)]

chart_style = st.selectbox("Chart Style", ["Table","North Indian","South Indian"], index=2)
render_svg = st.checkbox("Render charts as SVG (crispest)", value=False,
                         help="Vector rendering keeps lines ultra-thin at any zoom")

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

# -------- RENDER HELPERS (do NOT wrap in st.write) --------
def show_png(fig):
    fig.tight_layout(pad=0.12)
    st.pyplot(fig, use_container_width=False, dpi=300)

def show_svg(fig, width_px=240):
    buf = io.BytesIO()
    fig.savefig(buf, format="svg", bbox_inches="tight", pad_inches=0.02)
    st.image(buf.getvalue(), width=width_px)

# -------- OUTPUT --------
if st.session_state.chart_data:
    chart_data = st.session_state.chart_data
    st.markdown("---")
    st.markdown(f"""
    <div class="summary-box">
        <h3>Chart Summary</h3>
        <div class="summary-item"><strong>Name:</strong> {chart_data['name']}</div>
        <div class="summary-item"><strong>Lagna:</strong> {chart_data['lagna_sign']} ({chart_data['lagna_sid']:.2f}°)</div>
        <div class="summary-item"><strong>Rasi (Moon Sign):</strong> {chart_data['moon_rasi']}</div>
        <div class="summary-item"><strong>Nakshatra:</strong> {chart_data['moon_nakshatra']} (Pada {chart_data['moon_pada']})</div>
    </div>
    """, unsafe_allow_html=True)

    st.subheader("Planetary Positions")
    st.dataframe(chart_data['df_planets'], hide_index=True, use_container_width=True)

    # Side-by-side charts
    st.subheader("Rasi (D1) & Navamsa (D9)")
    if chart_style == "Table":
        ca, cb = st.columns(2, gap="small")
        with ca:
            st.markdown("**Rasi (D1)**")
            st.dataframe(chart_data['df_rasi'], hide_index=True, use_container_width=True)
        with cb:
            st.markdown("**Navamsa (D9)**")
            st.dataframe(chart_data['df_nav'], hide_index=True, use_container_width=True)
    else:
        # Build house->sign maps
        house_to_sign_rasi = {h: get_sign((chart_data['lagna_sid'] + (h-1)*30) % 360) for h in range(1,13)}
        house_to_sign_nav  = {h: get_sign((chart_data['nav_lagna'] + (h-1)*30) % 360) for h in range(1,13)}

        col1, col2 = st.columns(2, gap="small")
        size = (1.8, 1.8)  # smaller logical size

        with col1:
            if chart_style == "North Indian":
                fig, ax = plt.subplots(figsize=size)
                plot_north_indian_style(ax, chart_data['house_to_planets_rasi'], house_to_sign_rasi, 'Rasi Chart (North Indian)')
            else:
                fig, ax = plt.subplots(figsize=size)
                plot_south_indian_style(ax, chart_data['house_to_planets_rasi'], chart_data['lagna_sign'], 'Rasi Chart (South Indian)')
            show_svg(fig) if render_svg else show_png(fig)

        with col2:
            if chart_style == "North Indian":
                fig, ax = plt.subplots(figsize=size)
                plot_north_indian_style(ax, chart_data['house_to_planets_nav'], house_to_sign_nav, 'Navamsa Chart (North Indian)')
            else:
                fig, ax = plt.subplots(figsize=size)
                plot_south_indian_style(ax, chart_data['house_to_planets_nav'], chart_data['nav_lagna_sign'], 'Navamsa Chart (South Indian)')
            show_svg(fig) if render_svg else show_png(fig)

    st.subheader("House Analysis")
    st.dataframe(chart_data['df_house_status'], hide_index=True, use_container_width=True)

    st.subheader(f"Vimshottari Dasa ({chart_data['selected_depth']})")
    dasa_rows = []
    for lord, start, end, _subs in chart_data['dasa_periods_filtered']:
        dur = end - start
        dasa_rows.append({'Planet': lord, 'Start': start.strftime('%Y-%m-%d'),
                          'End': end.strftime('%Y-%m-%d'), 'Duration': duration_str(dur, 'dasa')})
    st.dataframe(pd.DataFrame(dasa_rows), hide_index=True, use_container_width=True)

    # Nested expanders (unchanged logic)
    max_depth = chart_data['max_depth']
    if max_depth >= 2:
        with st.expander("View Bhuktis (Sub-periods)", expanded=False):
            dp = chart_data['dasa_periods_filtered']
            if dp:
                dasa_options = [f"{p[0]} ({p[1].strftime('%Y-%m-%d')} - {p[2].strftime('%Y-%m-%d')})" for p in dp]
                selected_dasa = st.selectbox("Select Dasa:", dasa_options, key="dasa_select")
                sel_idx = dasa_options.index(selected_dasa)
                bhuktis = dp[sel_idx][3]
                bh_rows = []
                for b_lord, b_start, b_end, _ in bhuktis:
                    bh_rows.append({'Planet': b_lord, 'Start': b_start.strftime('%Y-%m-%d'),
                                    'End': b_end.strftime('%Y-%m-%d'),
                                    'Duration': duration_str(b_end-b_start, 'bhukti')})
                st.dataframe(pd.DataFrame(bh_rows), hide_index=True, use_container_width=True)

                if max_depth >= 3:
                    with st.expander("View Antharas", expanded=False):
                        if bhuktis:
                            bh_opts = [f"{p[0]} ({p[1].strftime('%Y-%m-%d')} - {p[2].strftime('%Y-%m-%d')})" for p in bhuktis]
                            selected_bh = st.selectbox("Select Bhukti:", bh_opts, key="bhukti_select")
                            sidx = bh_opts.index(selected_bh)
                            antharas = bhuktis[sidx][3]
                            an_rows = []
                            for a_lord, a_start, a_end, _ in antharas:
                                an_rows.append({'Planet': a_lord, 'Start': a_start.strftime('%Y-%m-%d %H:%M'),
                                                'End': a_end.strftime('%Y-%m-%d %H:%M'),
                                                'Duration': duration_str(a_end-a_start, 'anthara')})
                            st.dataframe(pd.DataFrame(an_rows), hide_index=True, use_container_width=True)

                            if max_depth >= 4:
                                with st.expander("View Sukshmas", expanded=False):
                                    if antharas:
                                        an_opts = [f"{p[0]} ({p[1].strftime('%Y-%m-%d %H:%M')} - {p[2].strftime('%Y-%m-%d %H:%M')})" for p in antharas]
                                        selected_an = st.selectbox("Select Anthara:", an_opts, key="anthara_select")
                                        aidx = an_opts.index(selected_an)
                                        sukshmas = antharas[aidx][3]
                                        sk_rows = []
                                        for s_lord, s_start, s_end, _ in sukshmas:
                                            sk_rows.append({'Planet': s_lord, 'Start': s_start.strftime('%Y-%m-%d %H:%M'),
                                                            'End': s_end.strftime('%Y-%m-%d %H:%M'),
                                                            'Duration': duration_str(s_end-s_start, 'sukshma')})
                                        st.dataframe(pd.DataFrame(sk_rows), hide_index=True, use_container_width=True)

                                        if max_depth >= 5:
                                            with st.expander("View Pranas", expanded=False):
                                                if sukshmas:
                                                    sk_opts = [f"{p[0]} ({p[1].strftime('%Y-%m-%d %H:%M')} - {p[2].strftime('%Y-%m-%d %H:%M')})" for p in sukshmas]
                                                    selected_sk = st.selectbox("Select Sukshma:", sk_opts, key="sukshma_select")
                                                    sidx2 = sk_opts.index(selected_sk)
                                                    pranas = sukshmas[sidx2][3]
                                                    pr_rows = []
                                                    for pr_lord, pr_start, pr_end, _ in pranas:
                                                        pr_rows.append({'Planet': pr_lord, 'Start': pr_start.strftime('%Y-%m-%d %H:%M'),
                                                                        'End': pr_end.strftime('%Y-%m-%d %H:%M'),
                                                                        'Duration': duration_str(pr_end-pr_start, 'prana')})
                                                    st.dataframe(pd.DataFrame(pr_rows), hide_index=True, use_container_width=True)

                                                    if max_depth >= 6:
                                                        with st.expander("View Sub-Pranas", expanded=False):
                                                            if pranas:
                                                                pr_opts = [f"{p[0]} ({p[1].strftime('%Y-%m-%d %H:%M')} - {p[2].strftime('%Y-%m-%d %H:%M')})" for p in pranas]
                                                                selected_pr = st.selectbox("Select Prana:", pr_opts, key="prana_select")
                                                                pidx = pr_opts.index(selected_pr)
                                                                sub_pranas = pranas[pidx][3]
                                                                sp_rows = []
                                                                for sp_lord, sp_start, sp_end, _ in sub_pranas:
                                                                    sp_rows.append({'Planet': sp_lord, 'Start': sp_start.strftime('%Y-%m-%d %H:%M'),
                                                                                    'End': sp_end.strftime('%Y-%m-%d %H:%M'),
                                                                                    'Duration': duration_str(sp_end-sp_start, 'sub_prana')})
                                                                st.dataframe(pd.DataFrame(sp_rows), hide_index=True, use_container_width=True)

    st.info("Note: Periods filtered from birth. Durations approximate.")
else:
    st.info("Enter details above and click 'Generate Chart' to begin")

st.markdown("---")
st.caption("Sivapathy Astrology Data Generator")
