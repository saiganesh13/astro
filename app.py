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
lords = lords_full  # Use full for output, but short for table

nak_names = [
    'Ashwini', 'Bharani', 'Krittika', 'Rohini', 'Mrigashira', 'Ardra', 'Punarvasu', 'Pushya', 'Ashlesha',
    'Magha', 'Purva Phalguni', 'Uttara Phalguni', 'Hasta', 'Chitra', 'Swati', 'Vishakha', 'Anuradha',
    'Jyeshta', 'Mula', 'Purva Ashadha', 'Uttara Ashadha', 'Shravana', 'Dhanishta', 'Shatabhisha',
    'Purva Bhadrapada', 'Uttara Bhadrapada', 'Revati'
]

years = [7, 20, 6, 10, 7, 18, 16, 19, 17] * 3

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
    d = jd - 2451545.0  # Days from J2000
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

def duration_str(delta):
    total_days = delta.total_seconds() / 86400  # More precise if needed, but approx
    years = int(total_days / 365.25)
    rem_days = total_days % 365.25
    months = int(rem_days / 30.4375)  # Avg month
    days = int(rem_days % 30.4375)
    if years + months + days == 0:
        return "Less than 1 day"
    return f"{years}y {months}m {days}d"

def compute_chart(date_text, time_text, ampm, lat, lon, tz_offset, max_depth):
    date_obj = datetime.strptime(date_text, '%d:%m:%Y')
    time_obj = datetime.strptime(time_text + ' ' + ampm, '%I:%M %p')
    local_dt = date_obj.replace(hour=time_obj.hour, minute=time_obj.minute, second=0, microsecond=0)

    utc_dt = local_dt - timedelta(hours=tz_offset)
    t = Time(utc_dt)
    jd = t.jd
    year = utc_dt.year
    ayan = get_lahiri_ayanamsa(year)

    # Planetary positions (tropical)
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
        for name, body in planet_bodies.items():
            p = get_body(body, t)
            ecl = p.transform_to(GeocentricTrueEcliptic())
            lon_trop[name] = ecl.lon.deg

    # Rahu/Ketu (mean nodes)
    d = jd - 2451545.0
    T = d / 36525.0
    omega = (125.04452 - 1934.136261 * T + 0.0020708 * T**2 + T**3 / 450000) % 360
    lon_trop['rahu'] = omega
    lon_trop['ketu'] = (omega + 180) % 360

    # Sidereal longitudes
    lon_sid = {p: get_sidereal_lon(lon_trop[p], ayan) for p in lon_trop}

    # Lagna (tropical asc, then sidereal)
    asc_trop = get_ascendant(jd, lat, lon)
    lagna_sid = get_sidereal_lon(asc_trop, ayan)

    # Planetary details
    data = []
    # Asc
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

    # Rasi Chart
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

    # Navamsa
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

    # Dasa - FIXED: correct dasa_start_dt
    moon_lon = lon_sid['moon']
    dasa_start_idx, balance_years = generate_vimshottari_dasa(moon_lon)
    full_first = years[dasa_start_idx]
    passed_years = full_first - balance_years
    dasa_start_dt = utc_dt - timedelta(days=passed_years * 365.25)

    # Generate full 120 years from dasa start
    dasa_periods = generate_periods(dasa_start_dt, dasa_start_idx, 120, level='dasa', max_depth=max_depth)
    dasa_periods_filtered = filter_from_birth(dasa_periods, utc_dt)

    max_depth_options = {1: 'Dasa only', 2: 'Dasa + Bhukti', 3: 'Dasa + Bhukti + Anthara', 
                         4: 'Dasa + Bhukti + Anthara + Sukshma', 5: 'Dasa + Bhukti + Anthara + Sukshma + Prana',
                         6: 'Dasa + Bhukti + Anthara + Sukshma + Prana + Sub-Prana'}
    selected_depth = max_depth_options[max_depth]

    lagna_sign = get_sign(lagna_sid)
    nav_lagna_sign = get_sign(nav_lagna)

    return {
        'df_planets': df_planets,
        'df_rasi': df_rasi,
        'df_nav': df_nav,
        'dasa_periods_filtered': dasa_periods_filtered,
        'lagna_sid': lagna_sid,
        'nav_lagna': nav_lagna,
        'lagna_sign': lagna_sign,
        'nav_lagna_sign': nav_lagna_sign,
        'selected_depth': selected_depth,
        'utc_dt': utc_dt,
        'max_depth': max_depth
    }

# Streamlit UI
st.set_page_config(page_title="Vedic Astrology", layout="wide")

# Custom CSS for styling - green only for button
st.markdown("""
<style>
    .stApp {
        background-color: white;
        color: #125336;
    }
    .stTextInput > div > div > input {
        background-color: white;
        color: #125336;
        border: 1px solid #125336;
    }
    .stTextInput > div > div > div > label {
        color: #125336;
    }
    .stSelectbox > label {
        color: #125336;
    }
    .stSelectbox > div > div > select {
        background-color: white;
        color: #125336;
        border: 1px solid #125336;
    }
    .stNumberInput > div > div > div > label {
        color: #125336;
    }
    .stNumberInput > div > div > input {
        background-color: white;
        color: #125336;
        border: 1px solid #125336;
    }
    .stButton > button {
        background-color: #125336;
        color: white;
        border: none;
    }
    .stButton > button:hover {
        background-color: #0a3d22;
    }
    .stMarkdown {
        color: #125336;
    }
    .stDataFrame {
        background-color: white;
        color: #125336;
    }
    h1, h2, h3 {
        color: #125336 !important;
    }
</style>
""", unsafe_allow_html=True)

st.title("🪐 Vedic Astrology Chart Generator")

# Initialize session state for chart data
if 'chart_data' not in st.session_state:
    st.session_state.chart_data = None

# Initialize geolocator with rate limiter
@st.cache_resource
def get_geolocator():
    geolocator = Nominatim(user_agent="vedic_astro_app")
    return RateLimiter(geolocator.geocode, min_delay_seconds=1)

geocode = get_geolocator()

# Inputs in columns
col1, col2 = st.columns(2)
with col1:
    date_text = st.text_input("Enter Date (DD:MM:YYYY):")
    time_text = st.text_input("Enter Time (HH:MM) in 12-hour format:")
    ampm = st.selectbox("AM/PM:", ["AM", "PM"])
with col2:
    # City search with suggestions
    city_query = st.text_input("Enter City Name (for lat/lon lookup):")
    location_data = None
    if city_query:
        try:
            # Search for multiple locations using geocode
            locations = geocode(city_query, exactly_one=False, limit=5)
            if locations:
                location_options = [f"{loc.address} (Lat: {loc.latitude:.2f}, Lon: {loc.longitude:.2f})" for loc in locations]
                selected_idx = st.selectbox("Select Location:", options=location_options, index=0, key="location_select")
                selected_location = locations[location_options.index(selected_idx)]
                location_data = {'lat': selected_location.latitude, 'lon': selected_location.longitude, 'address': selected_location.address}
                st.success(f"Selected: {selected_location.address}")
            else:
                st.warning("No locations found. Trying fallback...")
                city_key = city_query.title()
                if city_key in cities_fallback:
                    location_data = cities_fallback[city_key]
                    st.info("Using fallback data.")
        except Exception as e:
            st.error(f"Geocoding error: {e}")
            # Fallback
            city_key = city_query.title()
            if city_key in cities_fallback:
                location_data = cities_fallback[city_key]
                st.info("Using fallback data.")

    if location_data:
        lat = location_data['lat']
        lon = location_data['lon']
        st.write(f"Using location: Lat {lat:.2f}, Lon {lon:.2f}")
    else:
        lat = 13.08  # Default Chennai
        lon = 80.27
        st.warning("Using default location: Chennai")

    tz_offset = st.number_input("Timezone offset (hours, e.g., 5.5 for IST):", value=5.5, step=0.5)

col3, col4 = st.columns(2)
with col3:
    max_depth_options = {1: 'Dasa only', 2: 'Dasa + Bhukti', 3: 'Dasa + Bhukti + Anthara', 
                         4: 'Dasa + Bhukti + Anthara + Sukshma', 5: 'Dasa + Bhukti + Anthara + Sukshma + Prana',
                         6: 'Dasa + Bhukti + Anthara + Sukshma + Prana + Sub-Prana'}
    selected_depth_str = st.selectbox("Select max depth for periods:", options=list(max_depth_options.values()), index=2)
    max_depth = list(max_depth_options.keys())[list(max_depth_options.values()).index(selected_depth_str)]

# Button to generate or regenerate
if st.button("Generate Chart", use_container_width=True):
    st.session_state.chart_data = compute_chart(date_text, time_text, ampm, lat, lon, tz_offset, max_depth)
    st.rerun()

# Display chart if computed
if st.session_state.chart_data:
    chart_data = st.session_state.chart_data
    st.subheader("=== PLANETARY DETAILS ===")
    st.table(chart_data['df_planets'])

    st.subheader("=== RASI CHART ===")
    st.write(f"Lagna: {chart_data['lagna_sign']} ({chart_data['lagna_sid']:.2f}°)")
    st.table(chart_data['df_rasi'])

    st.subheader("=== NAVAMSA CHART ===")
    st.write(f"Navamsa Lagna: {chart_data['nav_lagna_sign']} ({chart_data['nav_lagna']:.2f}°)")
    st.table(chart_data['df_nav'])

    st.subheader(f"=== VIMSHOTTARI {chart_data['selected_depth'].upper()} ===")

    dasa_periods_filtered = chart_data['dasa_periods_filtered']
    utc_dt = chart_data['utc_dt']
    max_depth = chart_data['max_depth']

    # Dasa Table
    dasa_data = []
    for lord, start, end, _ in dasa_periods_filtered:
        dur = end - start
        dasa_data.append({'Planet': lord, 'Start': start.strftime('%Y-%m-%d'), 'End': end.strftime('%Y-%m-%d'), 'Duration': duration_str(dur)})
    df_dasa = pd.DataFrame(dasa_data)
    st.table(df_dasa)

    # Nested for deeper levels using expanders and selects
    if max_depth >= 2:
        with st.expander("Bhuktis (Select Dasa to view)"):
            if dasa_periods_filtered:
                dasa_options = {i: f"{p[0]} ({p[1].strftime('%Y-%m-%d')} - {p[2].strftime('%Y-%m-%d')})" for i, p in enumerate(dasa_periods_filtered)}
                selected_dasa_idx = st.selectbox("Select Dasa:", options=list(dasa_options.values()), index=0, format_func=lambda x: x)
                sel_idx = list(dasa_options.keys())[list(dasa_options.values()).index(selected_dasa_idx)]
                selected_dasa_period = dasa_periods_filtered[sel_idx]
                bhuktis = selected_dasa_period[3]
                bhukti_data = []
                for b_lord, b_start, b_end, _ in bhuktis:
                    dur = b_end - b_start
                    bhukti_data.append({'Planet': b_lord, 'Start': b_start.strftime('%Y-%m-%d'), 'End': b_end.strftime('%Y-%m-%d'), 'Duration': duration_str(dur)})
                df_bhukti = pd.DataFrame(bhukti_data)
                st.table(df_bhukti)

                if max_depth >= 3:
                    with st.expander("Antharas (Select Bhukti to view)"):
                        if bhuktis:
                            bhukti_options = {j: f"{p[0]} ({p[1].strftime('%Y-%m-%d')} - {p[2].strftime('%Y-%m-%d')})" for j, p in enumerate(bhuktis)}
                            selected_bhukti_idx = st.selectbox("Select Bhukti:", options=list(bhukti_options.values()), index=0, format_func=lambda x: x)
                            sel_b_idx = list(bhukti_options.keys())[list(bhukti_options.values()).index(selected_bhukti_idx)]
                            selected_bhukti_period = bhuktis[sel_b_idx]
                            antharas = selected_bhukti_period[3]
                            anthara_data = []
                            for a_lord, a_start, a_end, _ in antharas:
                                dur = a_end - a_start
                                anthara_data.append({'Planet': a_lord, 'Start': a_start.strftime('%Y-%m-%d %H:%M'), 'End': a_end.strftime('%Y-%m-%d %H:%M'), 'Duration': duration_str(dur)})
                            df_anthara = pd.DataFrame(anthara_data)
                            st.table(df_anthara)

                            if max_depth >= 4:
                                with st.expander("Sukshmas (Select Anthara to view)"):
                                    if antharas:
                                        anthara_options = {k: f"{p[0]} ({p[1].strftime('%Y-%m-%d %H:%M')} - {p[2].strftime('%Y-%m-%d %H:%M')})" for k, p in enumerate(antharas)}
                                        selected_anthara_idx = st.selectbox("Select Anthara:", options=list(anthara_options.values()), index=0, format_func=lambda x: x)
                                        sel_a_idx = list(anthara_options.keys())[list(anthara_options.values()).index(selected_anthara_idx)]
                                        selected_anthara_period = antharas[sel_a_idx]
                                        sukshmas = selected_anthara_period[3]
                                        sukshma_data = []
                                        for s_lord, s_start, s_end, _ in sukshmas:
                                            dur = s_end - s_start
                                            sukshma_data.append({'Planet': s_lord, 'Start': s_start.strftime('%Y-%m-%d %H:%M'), 'End': s_end.strftime('%Y-%m-%d %H:%M'), 'Duration': duration_str(dur)})
                                        df_sukshma = pd.DataFrame(sukshma_data)
                                        st.table(df_sukshma)

                                        if max_depth >= 5:
                                            with st.expander("Pranas (Select Sukshma to view)"):
                                                if sukshmas:
                                                    sukshma_options = {l: f"{p[0]} ({p[1].strftime('%Y-%m-%d %H:%M')} - {p[2].strftime('%Y-%m-%d %H:%M')})" for l, p in enumerate(sukshmas)}
                                                    selected_sukshma_idx = st.selectbox("Select Sukshma:", options=list(sukshma_options.values()), index=0, format_func=lambda x: x)
                                                    sel_s_idx = list(sukshma_options.keys())[list(sukshma_options.values()).index(selected_sukshma_idx)]
                                                    selected_sukshma_period = sukshmas[sel_s_idx]
                                                    pranas = selected_sukshma_period[3]
                                                    prana_data = []
                                                    for pr_lord, pr_start, pr_end, _ in pranas:
                                                        dur = pr_end - pr_start
                                                        prana_data.append({'Planet': pr_lord, 'Start': pr_start.strftime('%Y-%m-%d %H:%M'), 'End': pr_end.strftime('%Y-%m-%d %H:%M'), 'Duration': duration_str(dur)})
                                                    df_prana = pd.DataFrame(prana_data)
                                                    st.table(df_prana)

                                                    if max_depth >= 6:
                                                        with st.expander("Sub-Pranas (Select Prana to view)"):
                                                            if pranas:
                                                                prana_options = {m: f"{p[0]} ({p[1].strftime('%Y-%m-%d %H:%M')} - {p[2].strftime('%Y-%m-%d %H:%M')})" for m, p in enumerate(pranas)}
                                                                selected_prana_idx = st.selectbox("Select Prana:", options=list(prana_options.values()), index=0, format_func=lambda x: x)
                                                                sel_pr_idx = list(prana_options.keys())[list(prana_options.values()).index(selected_prana_idx)]
                                                                selected_prana_period = pranas[sel_pr_idx]
                                                                sub_pranas = selected_prana_period[3]
                                                                sub_prana_data = []
                                                                for sp_lord, sp_start, sp_end, _ in sub_pranas:
                                                                    dur = sp_end - sp_start
                                                                    sub_prana_data.append({'Planet': sp_lord, 'Start': sp_start.strftime('%Y-%m-%d %H:%M'), 'End': sp_end.strftime('%Y-%m-%d %H:%M'), 'Duration': duration_str(dur)})
                                                                df_sub_prana = pd.DataFrame(sub_prana_data)
                                                                st.table(df_sub_prana)
    st.info("Note: Periods filtered from birth. Durations approximate. Deeper levels use nested expanders for navigation.")
else:
    st.info("Enter details and click 'Generate Chart' to begin.")

# To host: Save as app.py, run `streamlit run app.py` locally. For cloud hosting, push to GitHub and deploy on Streamlit Cloud.
