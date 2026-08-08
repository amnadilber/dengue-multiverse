"""
Representative climate locations for each reporting unit.

Every window needs one point at which to sample temperature and rainfall. There
is no principled single choice — a country's climate is not a point — so the rule
applied throughout is the same one used for Pakistan: **the largest reporting
centre**, which for a national series is normally the capital or the largest
city, and for a subnational unit its administrative seat.

This is a simplification with consequences, and they are measured rather than
assumed: the Pakistani sensitivity analysis found that swapping Lahore for
Peshawar moved the national R0 estimate by 0.32, more than the bootstrap
sampling interval. The same test is applied across the global sample.

Coordinates are decimal degrees, from standard gazetteer values for each city.
They are recorded here as data rather than fetched from a geocoding service so
that the analysis has no network dependency at run time and the exact points used
are auditable.
"""

from __future__ import annotations

# --- National reporting units: capital or largest city ---------------------
COUNTRY_POINTS: dict[str, tuple[str, float, float]] = {
    "ARGENTINA": ("Buenos Aires", -34.61, -58.38),
    "BARBADOS": ("Bridgetown", 13.10, -59.62),
    "BELIZE": ("Belize City", 17.50, -88.20),
    "BOLIVIA": ("Santa Cruz de la Sierra", -17.78, -63.18),
    "BRAZIL": ("Sao Paulo", -23.55, -46.63),
    "CAMBODIA": ("Phnom Penh", 11.56, 104.92),
    "COLOMBIA": ("Bogota", 4.71, -74.07),
    "COSTA RICA": ("San Jose", 9.93, -84.08),
    "DOMINICAN REPUBLIC": ("Santo Domingo", 18.49, -69.93),
    "ECUADOR": ("Guayaquil", -2.19, -79.89),
    "EL SALVADOR": ("San Salvador", 13.69, -89.19),
    "FRENCH GUIANA": ("Cayenne", 4.93, -52.33),
    "GRENADA": ("Saint George's", 12.06, -61.75),
    "GUADELOUPE": ("Pointe-a-Pitre", 16.24, -61.53),
    "GUATEMALA": ("Guatemala City", 14.63, -90.51),
    "GUYANA": ("Georgetown", 6.80, -58.16),
    "HONDURAS": ("Tegucigalpa", 14.07, -87.19),
    "JAMAICA": ("Kingston", 17.98, -76.79),
    "LAO PEOPLE'S DEMOCRATIC REPUBLIC": ("Vientiane", 17.97, 102.63),
    "MEXICO": ("Mexico City", 19.43, -99.13),
    "NEW CALEDONIA": ("Noumea", -22.28, 166.46),
    "NICARAGUA": ("Managua", 12.11, -86.24),
    "PAKISTAN": ("Lahore", 31.55, 74.35),
    "PANAMA": ("Panama City", 8.98, -79.52),
    "PARAGUAY": ("Asuncion", -25.26, -57.58),
    "PERU": ("Lima", -12.05, -77.04),
    "PHILIPPINES": ("Manila", 14.60, 120.98),
    "PUERTO RICO": ("San Juan", 18.47, -66.11),
    "SAINT BARTHELEMY": ("Gustavia", 17.90, -62.85),
    "SAINT MARTIN": ("Marigot", 18.07, -63.08),
    "SINGAPORE": ("Singapore", 1.35, 103.82),
    "UNITED STATES OF AMERICA": ("San Juan PR / Key West", 24.56, -81.78),
    "VENEZUELA": ("Caracas", 10.49, -66.88),
    "VIET NAM": ("Ho Chi Minh City", 10.82, 106.63),
    "YEMEN": ("Aden", 12.79, 45.04),
}

# --- Subnational units: administrative seat --------------------------------
# Only the units that actually produce a usable window are listed. Names are the
# uppercased forms OpenDengue uses, which do not always match the official
# spelling.
ADMIN1_POINTS: dict[tuple[str, str], tuple[str, float, float]] = {
    # Bolivia — departments
    ("BOLIVIA", "SANTA CRUZ"): ("Santa Cruz de la Sierra", -17.78, -63.18),
    ("BOLIVIA", "BENI"): ("Trinidad", -14.83, -64.90),
    ("BOLIVIA", "COCHABAMBA"): ("Cochabamba", -17.39, -66.16),
    ("BOLIVIA", "LA PAZ"): ("La Paz", -16.50, -68.15),
    ("BOLIVIA", "TARIJA"): ("Tarija", -21.53, -64.73),
    ("BOLIVIA", "CHUQUISACA"): ("Sucre", -19.03, -65.26),
    ("BOLIVIA", "PANDO"): ("Cobija", -11.03, -68.77),
    # Nicaragua — departments and autonomous regions
    ("NICARAGUA", "MANAGUA"): ("Managua", 12.11, -86.24),
    ("NICARAGUA", "CHINANDEGA"): ("Chinandega", 12.63, -87.13),
    ("NICARAGUA", "LEON"): ("Leon", 12.44, -86.88),
    ("NICARAGUA", "MASAYA"): ("Masaya", 11.97, -86.09),
    ("NICARAGUA", "GRANADA"): ("Granada", 11.93, -85.96),
    ("NICARAGUA", "CARAZO"): ("Jinotepe", 11.85, -86.20),
    ("NICARAGUA", "RIVAS"): ("Rivas", 11.44, -85.83),
    ("NICARAGUA", "CHONTALES"): ("Juigalpa", 12.11, -85.36),
    ("NICARAGUA", "MATAGALPA"): ("Matagalpa", 12.93, -85.92),
    ("NICARAGUA", "JINOTEGA"): ("Jinotega", 13.09, -86.00),
    ("NICARAGUA", "ESTELI"): ("Esteli", 13.09, -86.35),
    ("NICARAGUA", "MADRIZ"): ("Somoto", 13.48, -86.58),
    ("NICARAGUA", "NUEVA SEGOVIA"): ("Ocotal", 13.63, -86.48),
    ("NICARAGUA", "BOACO"): ("Boaco", 12.47, -85.66),
    ("NICARAGUA", "RIO SAN JUAN"): ("San Carlos", 11.12, -84.78),
    ("NICARAGUA", "RAAN"): ("Puerto Cabezas", 14.04, -83.39),
    ("NICARAGUA", "RAAS"): ("Bluefields", 12.01, -83.76),
    ("NICARAGUA", "ATLANTICO NORTE"): ("Puerto Cabezas", 14.04, -83.39),
    ("NICARAGUA", "ATLANTICO SUR"): ("Bluefields", 12.01, -83.76),
    # Mexico — states
    ("MEXICO", "VERACRUZ"): ("Veracruz", 19.17, -96.13),
    ("MEXICO", "GUERRERO"): ("Acapulco", 16.86, -99.88),
    ("MEXICO", "OAXACA"): ("Oaxaca", 17.07, -96.72),
    ("MEXICO", "CHIAPAS"): ("Tuxtla Gutierrez", 16.75, -93.12),
    ("MEXICO", "JALISCO"): ("Guadalajara", 20.67, -103.35),
    ("MEXICO", "MICHOACAN"): ("Morelia", 19.70, -101.19),
    ("MEXICO", "SINALOA"): ("Culiacan", 24.81, -107.39),
    ("MEXICO", "SONORA"): ("Hermosillo", 29.07, -110.96),
    ("MEXICO", "TAMAULIPAS"): ("Tampico", 22.25, -97.87),
    ("MEXICO", "YUCATAN"): ("Merida", 20.97, -89.62),
    ("MEXICO", "QUINTANA ROO"): ("Chetumal", 18.50, -88.30),
    ("MEXICO", "COLIMA"): ("Colima", 19.24, -103.72),
    ("MEXICO", "NAYARIT"): ("Tepic", 21.51, -104.89),
    ("MEXICO", "MORELOS"): ("Cuernavaca", 18.92, -99.23),
    ("MEXICO", "TABASCO"): ("Villahermosa", 17.99, -92.93),
    ("MEXICO", "SAN LUIS POTOSI"): ("San Luis Potosi", 22.16, -100.98),
    ("MEXICO", "PUEBLA"): ("Puebla", 19.04, -98.20),
    ("MEXICO", "NUEVO LEON"): ("Monterrey", 25.69, -100.32),
    # Peru — regions
    ("PERU", "LIMA"): ("Lima", -12.05, -77.04),
    ("PERU", "PIURA"): ("Piura", -5.19, -80.63),
    ("PERU", "LAMBAYEQUE"): ("Chiclayo", -6.77, -79.84),
    ("PERU", "LA LIBERTAD"): ("Trujillo", -8.11, -79.03),
    ("PERU", "LORETO"): ("Iquitos", -3.75, -73.25),
    ("PERU", "UCAYALI"): ("Pucallpa", -8.38, -74.55),
    ("PERU", "SAN MARTIN"): ("Tarapoto", -6.49, -76.37),
    ("PERU", "JUNIN"): ("Huancayo", -12.07, -75.21),
    ("PERU", "MADRE DE DIOS"): ("Puerto Maldonado", -12.60, -69.19),
    ("PERU", "CAJAMARCA"): ("Cajamarca", -7.16, -78.51),
    ("PERU", "TUMBES"): ("Tumbes", -3.57, -80.46),
    ("PERU", "ANCASH"): ("Chimbote", -9.08, -78.59),
    ("PERU", "ICA"): ("Ica", -14.07, -75.73),
    # Colombia — departments
    ("COLOMBIA", "ANTIOQUIA"): ("Medellin", 6.24, -75.58),
    ("COLOMBIA", "VALLE DEL CAUCA"): ("Cali", 3.45, -76.53),
    ("COLOMBIA", "SANTANDER"): ("Bucaramanga", 7.12, -73.13),
    ("COLOMBIA", "NORTE DE SANTANDER"): ("Cucuta", 7.89, -72.50),
    ("COLOMBIA", "TOLIMA"): ("Ibague", 4.44, -75.24),
    ("COLOMBIA", "HUILA"): ("Neiva", 2.93, -75.28),
    ("COLOMBIA", "META"): ("Villavicencio", 4.15, -73.63),
    ("COLOMBIA", "CUNDINAMARCA"): ("Girardot", 4.30, -74.80),
    ("COLOMBIA", "CESAR"): ("Valledupar", 10.46, -73.25),
    ("COLOMBIA", "CASANARE"): ("Yopal", 5.34, -72.40),
    # Dominican Republic and Ecuador
    ("DOMINICAN REPUBLIC", "DISTRITO NACIONAL"): ("Santo Domingo", 18.49, -69.93),
    ("DOMINICAN REPUBLIC", "SANTIAGO"): ("Santiago", 19.45, -70.70),
    ("ECUADOR", "GUAYAS"): ("Guayaquil", -2.19, -79.89),
    ("ECUADOR", "MANABI"): ("Portoviejo", -1.05, -80.45),
    ("ECUADOR", "LOS RIOS"): ("Babahoyo", -1.80, -79.53),
    ("ECUADOR", "EL ORO"): ("Machala", -3.26, -79.96),
    ("COLOMBIA", "ARAUCA"): ("Arauca", 7.09, -70.76),
    ("COLOMBIA", "CAQUETA"): ("Florencia", 1.61, -75.61),
    ("COLOMBIA", "CORDOBA"): ("Monteria", 8.75, -75.88),
    ("COLOMBIA", "CUNDIMARCA"): ("Girardot", 4.30, -74.80),
    ("COLOMBIA", "QUINDIO"): ("Armenia", 4.53, -75.68),
    ("COLOMBIA", "SUCRE"): ("Sincelejo", 9.30, -75.40),
    ("DOMINICAN REPUBLIC", "SAN CRISTOBAL"): ("San Cristobal", 18.42, -70.11),
    ("ECUADOR", "ORELLANA"): ("Coca", -0.47, -76.99),
    ("MEXICO", "BAJA CALIFORNIA SUR"): ("La Paz", 24.14, -110.31),
    ("NICARAGUA", "REGION AUTONOMA DEL ATLANTICO SUR"): ("Bluefields", 12.01, -83.76),
    ("PERU", "AMAZONAS"): ("Chachapoyas", -6.23, -77.87),
}


def point_for(country: str, unit: str, level: str):
    """Return (name, lat, lon) for a reporting unit, or None if unknown."""
    country = country.upper().strip()
    if level == "national":
        return COUNTRY_POINTS.get(country)
    return ADMIN1_POINTS.get((country, unit.upper().strip()))


def coverage(windows) -> tuple[int, int, list]:
    """How many windows can be given a climate point? Returns (have, total, missing)."""
    missing = []
    have = 0
    for w in windows:
        if point_for(w["country"], w["unit"], w["level"]):
            have += 1
        else:
            missing.append((w["country"], w["unit"], w["level"]))
    return have, len(windows), sorted(set(missing))


#: Latitude offset, in degrees, defining the alternative climate point used to
#: test whether the choice of grid cell changes a conclusion. One degree is
#: about 110 km — deliberately modest, since it is well inside the area any
#: national or provincial case series aggregates over, so a verdict that changes
#: at this distance is changing for a reason no analyst could defend.
CLIMATE_OFFSET_DEG = 1.0


def offset_point(lat: float, lon: float) -> tuple[float, float]:
    """The alternative climate point for a location.

    Shifts toward the equator, which at every latitude represented in this study
    keeps the result on the globe and avoids having to reflect the offset at a
    pole. Lives here rather than in a pipeline script because two steps need it
    and they must agree exactly: a mismatch would silently look for a file that
    was never downloaded and report the location as unavailable.
    """
    shifted = lat - CLIMATE_OFFSET_DEG if lat > 0 else lat + CLIMATE_OFFSET_DEG
    return round(shifted, 4), round(lon, 4)
