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
import io

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
    sub  = (star + int((pos/dnak)*9)) % 9
    return nak_names[idx], pada, lords_short[star], lords_short[sub]

def compute_sidereal_positions(utc_dt):
    t = Time(utc_dt); jd = t.jd; ayan = get_lahiri_ayanamsa(utc_dt.year)

    with solar_system_ephemeris.set('builtin'):
        lon_trop = {}
        for nm in ['sun','moon','mercury','venus','mars','jupiter','saturn']:
            ecl = get_body(nm, t).transform_to(GeocentricTrueEcliptic()); lon_trop[nm] = ecl.lon.deg
    d = jd - 2451545.0; T = d/36525.0
    omega = (125.04452 - 1934.136261*T + 0.0020708*T**2 + T**3/450000) % 360
    lon_trop['rahu'] = omega; lon_trop['ketu'] = (omega + 180) % 360

    lon_sid = {p: get_sidereal_lon(lon_trop[p], ayan) for p in lon_trop}
    return lon_sid

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

def compute_chart(name, date_obj, time_str, lat, lon, tz_offset, max_depth):
    # parse time
    try:
        hour, minute = map(int, time_str.split(':'))
        if not (0<=hour<=23 and 0<=minute<=59): raise ValueError
    except:
        raise ValueError("Time must be in HH:MM format (24-hour)")
    local_dt = datetime.combine(date_obj, datetime.min.time().replace(hour=hour, minute=minute))
    utc_dt = local_dt - timedelta(hours=tz_offset)
    lon_sid = compute_sidereal_positions(utc_dt)
    jd = Time(utc_dt).jd
    ayan = get_lahiri_ayanamsa(utc_dt.year)
    lagna_sid = get_sidereal_lon(get_ascendant(jd, lat, lon), ayan)

    # planets table
    rows = []
    asc_deg = lagna_sid % 360; asc_sign = get_sign(asc_deg)
    a_nak, a_pada, a_ld, a_sl = get_nakshatra_details(asc_deg)
    rows.append(['Asc', f"{asc_deg:.2f}", asc_sign, a_nak, a_pada, f"{a_ld}/{a_sl}"])
    for p in ['sun','moon','mars','mercury','jupiter','venus','saturn','rahu','ketu']:
        L = lon_sid[p]; sign = get_sign(L); nak, pada, ld, sl = get_nakshatra_details(L)
        rows.append([p.capitalize(), f"{L:.2f}", sign, nak, pada, f"{ld}/{sl}"])
    df_planets = pd.DataFrame(rows, columns=['Planet','Deg','Sign','Nakshatra','Pada','Ld/SL'])

    # rasi houses
    house_planets_rasi = defaultdict(list); positions = {**lon_sid,'asc':lagna_sid}
    for p,L in positions.items():
        house_planets_rasi[get_house(L, lagna_sid)].append(p.capitalize() if p!='asc' else 'Asc')
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

    # aspects table (unchanged logic)
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
        'house_to_planets_rasi': house_planets_rasi, 'house_to_planets_nav': house_planets_nav,
        'natal_moon_lon': moon_lon, 'tz_offset': tz_offset
    }

# ---- South Indian plotter ONLY (extra top padding) ----
def plot_south_indian_style(ax, house_to_planets, lagna_sign, title):
    sign_positions = {'Pisces':(0,3),'Aries':(1,3),'Taurus':(2,3),'Gemini':(3,3),
                      'Cancer':(3,2),'Leo':(3,1),'Virgo':(3,0),
                      'Libra':(2,0),'Scorpio':(1,0),'Sagittarius':(0,0),
                      'Capricorn':(0,1),'Aquarius':(0,2)}
    lagna_idx = sign_names.index(lagna_sign)
    house_for_sign = {s: ((i - lagna_idx) % 12) + 1 for i, s in enumerate(sign_names)}

    # tuned layout
    box_w, box_h, spacing, pad = 0.46, 0.46, 0.52, 0.02
    top_pad_extra = 0.020   # EXTRA space above the first planet line (increased)
    line_h_min    = 0.042   # taller minimum step for breathing room
    line_h_max    = 0.058
    planet_font   = 2.45    # small to support larger spacing

    for sign,(gx,gy) in sign_positions.items():
        h = house_for_sign[sign]
        planets = sorted(house_to_planets.get(h,[]))
        x = gx*spacing + 0.22; y = (3-gy)*spacing + 0.22

        ax.add_patch(FancyBboxPatch((x,y), box_w, box_h,
                                    boxstyle="round,pad=0.004",
                                    ec="black", fc="#F5F5F5",
                                    alpha=0.92, linewidth=0.32))
        # tiny sign label left-top
        ax.text(x+pad, y+pad, sign[:3], ha='left', va='top', fontsize=2.7)

        if planets:
            # reduce usable height to force looser stacking; add extra top margin
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

# Inputs
st.subheader("Birth Details")
name = st.text_input("Name", placeholder="Enter full name")
c1, c2, c3 = st.columns(3)
with c1:
    birth_date = st.date_input("Birth Date", value=datetime.now().date(),
                               min_value=datetime(1900,1,1).date(), max_value=datetime.now().date())
with c2:
    birth_time = st.text_input("Birth Time (HH:MM in 24-hour format)", placeholder="14:30")
with c3:
    tz_offset = st.number_input("Timezone offset (hrs)", value=5.5, step=0.5,
                                help="Offset from UTC (e.g., IST = 5.5)")

use_custom_coords = st.checkbox("Custom latitude and longitude?")
if use_custom_coords:
    clat, clon = st.columns(2)
    with clat: lat = st.number_input("Latitude", value=13.08, format="%.4f")
    with clon: lon = st.number_input("Longitude", value=80.27, format="%.4f")
else:
    city_query = st.text_input("Search City", placeholder="Start typing city name...", key="city_input")
    if city_query and len(city_query) >= 2:
        try:
            locations = geocode(city_query, exactly_one=False, limit=5)
            st.session_state.search_results = [{'display': loc.address, 'lat': loc.latitude, 'lon': loc.longitude, 'address': loc.address} for loc in (locations or [])]
        except:
            st.session_state.search_results = []
    else:
        st.session_state.search_results = []
    if st.session_state.search_results:
        opts = [r['display'] for r in st.session_state.search_results]
        sel = st.selectbox("Select location", options=opts, key="location_selector")
        i = opts.index(sel)
        lat = st.session_state.search_results[i]['lat']; lon = st.session_state.search_results[i]['lon']
        st.success(f"Selected: {st.session_state.search_results[i]['address']} (Lat: {lat:.2f}, Lon: {lon:.2f})")
    else:
        city_key = (city_query or "").title()
        if city_key in cities_fallback:
            lat = cities_fallback[city_key]['lat']; lon = cities_fallback[city_key]['lon']
        else:
            lat, lon = 13.08, 80.27
        st.info("Using default location: Chennai, India")

max_depth_options = {
    1:'Dasa only',2:'Dasa + Bhukti',3:'Dasa + Bhukti + Anthara',
    4:'Dasa + Bhukti + Anthara + Sukshma',5:'Dasa + Bhukti + Anthara + Sukshma + Prana',
    6:'Dasa + Bhukti + Anthara + Sukshma + Prana + Sub-Prana'
}
selected_depth_str = st.selectbox("Period Depth", list(max_depth_options.values()), index=2)
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

# render helpers
def show_svg(fig, width_px=240):
    buf = io.BytesIO()
    fig.savefig(buf, format="svg", bbox_inches="tight", pad_inches=0.02)
    st.image(buf.getvalue(), width=width_px)

def show_png(fig):
    fig.tight_layout(pad=0.10)
    st.pyplot(fig, use_container_width=False, dpi=300)

# Output
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

    # South Indian charts side-by-side
    st.subheader("Rasi (D1) & Navamsa (D9) — South Indian")
    col1, col2 = st.columns(2, gap="small")
    size = (1.8, 1.8)

    with col1:
        fig, ax = plt.subplots(figsize=size)
        plot_south_indian_style(ax, cd['house_to_planets_rasi'], cd['lagna_sign'], 'Rasi Chart (South Indian)')
        show_png(fig)

    with col2:
        fig, ax = plt.subplots(figsize=size)
        plot_south_indian_style(ax, cd['house_to_planets_nav'], cd['nav_lagna_sign'], 'Navamsa Chart (South Indian)')
        show_png(fig)

    st.subheader("House Analysis")
    st.dataframe(cd['df_house_status'], hide_index=True, use_container_width=True)

    # Today's Rasi Effects
    utc_now = datetime.utcnow()
    tz_offset = cd['tz_offset']
    local_now = utc_now + timedelta(hours=tz_offset)
    today_date_obj = local_now.date()
    today_time_str = local_now.strftime('%H:%M')
    hour, minute = map(int, today_time_str.split(':'))
    local_dt_today = datetime.combine(today_date_obj, datetime.min.time().replace(hour=hour, minute=minute))
    utc_dt_today = local_dt_today - timedelta(hours=tz_offset)
    lon_sid_today = compute_sidereal_positions(utc_dt_today)

    chandra_lagna_lon = cd['natal_moon_lon']
    transit_planet_to_house = {p.capitalize(): get_house(lon_sid_today[p], chandra_lagna_lon) for p in lon_sid_today}
    house_planets_transit = defaultdict(list)
    for p, h in transit_planet_to_house.items():
        house_planets_transit[h].append(p)
    aspects_dict = {'Sun':[7],'Moon':[7],'Mars':[4,7,8],'Mercury':[7],'Jupiter':[5,7,9],'Venus':[7],'Saturn':[3,7,10]}
    house_status_today = []
    for h in range(1,13):
        house_sign_deg = (chandra_lagna_lon + (h-1)*30 ) % 360
        house_sign = get_sign(house_sign_deg)
        lord = sign_lords[sign_names.index(house_sign)]
        planets = sorted(house_planets_transit[h])
        asp = []
        for planet, offs in aspects_dict.items():
            if planet in transit_planet_to_house:
                ph = transit_planet_to_house[planet]
                for off in offs:
                    if ((ph-1 + (off-1)) % 12 ) + 1 == h:
                        asp.append(planet)
        lord_house = get_house(lon_sid_today[lord.lower()], chandra_lagna_lon)
        house_status_today.append([f"House {h}",
                                   ', '.join(planets) if planets else 'Empty',
                                   ', '.join(sorted(asp)) if asp else 'None',
                                   lord,
                                   f"House {lord_house}"])
    df_transit_house = pd.DataFrame(house_status_today, columns=['House','Planets (Transit)','Aspects from','Lord','Lord in (Transit)'])

    st.subheader("Today's Rasi Effects (Chandra Lagna)")
    st.dataframe(df_transit_house, hide_index=True, use_container_width=True)
    st.caption(f"Based on transits as of {local_now.strftime('%Y-%m-%d %H:%M')} local time")

    # Vimshottari with full nested expanders
    st.subheader(f"Vimshottari Dasa ({cd['selected_depth']})")
    dasa_rows = [{'Planet': lord, 'Start': s.strftime('%Y-%m-%d'),
                  'End': e.strftime('%Y-%m-%d'), 'Duration': duration_str(e-s,'dasa')}
                 for lord, s, e, _ in cd['dasa_periods_filtered']]
    st.dataframe(pd.DataFrame(dasa_rows), hide_index=True, use_container_width=True)

    max_depth = cd['max_depth']
    dp = cd['dasa_periods_filtered']
    if max_depth >= 2:
        with st.expander("View Bhuktis (Sub-periods)", expanded=False):
            if dp:
                d_opt = [f"{p[0]} ({p[1].strftime('%Y-%m-%d')} - {p[2].strftime('%Y-%m-%d')})" for p in dp]
                sel = st.selectbox("Select Dasa:", d_opt, key="dasa_select")
                bhuktis = dp[d_opt.index(sel)][3]
                st.dataframe(pd.DataFrame(
                    [{'Planet': l, 'Start': s.strftime('%Y-%m-%d'), 'End': e.strftime('%Y-%m-%d'),
                      'Duration': duration_str(e-s,'bhukti')} for l,s,e,_ in bhuktis]),
                    hide_index=True, use_container_width=True)
                if max_depth >= 3:
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
                            if max_depth >= 4:
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
                                        if max_depth >= 5:
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
                                                    if max_depth >= 6:
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
else:
    st.info("Enter details above and click 'Generate Chart' to begin")

st.markdown("---")
st.caption("Sivapathy Astrology Data Generator")
