// Copyright (C) 2026 CRV Science Network (https://crvscience.com)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://gnu.org>.

/// 12x8 US State + DC Cartogram Grid Layout `(col, row, fips, abbr)` for zero-asset native GPU rendering.
pub const US_CARTOGRAM_CELLS: &[(usize, usize, &str, &str)] = &[
    (0, 0, "02", "AK"), (11, 0, "23", "ME"),
    (10, 1, "50", "VT"), (11, 1, "33", "NH"),
    (1, 2, "53", "WA"), (2, 2, "16", "ID"), (3, 2, "30", "MT"), (4, 2, "38", "ND"),
    (5, 2, "27", "MN"), (6, 2, "55", "WI"), (7, 2, "26", "MI"), (9, 2, "36", "NY"),
    (10, 2, "25", "MA"), (11, 2, "44", "RI"),
    (1, 3, "41", "OR"), (2, 3, "32", "NV"), (3, 3, "56", "WY"), (4, 3, "46", "SD"),
    (5, 3, "19", "IA"), (6, 3, "17", "IL"), (7, 3, "18", "IN"), (8, 3, "39", "OH"),
    (9, 3, "42", "PA"), (10, 3, "34", "NJ"), (11, 3, "09", "CT"),
    (1, 4, "06", "CA"), (2, 4, "49", "UT"), (3, 4, "08", "CO"), (4, 4, "31", "NE"),
    (5, 4, "29", "MO"), (6, 4, "21", "KY"), (7, 4, "54", "WV"), (8, 4, "51", "VA"),
    (9, 4, "24", "MD"), (10, 4, "10", "DE"),
    (2, 5, "04", "AZ"), (3, 5, "35", "NM"), (4, 5, "20", "KS"), (5, 5, "05", "AR"),
    (6, 5, "47", "TN"), (7, 5, "37", "NC"), (8, 5, "45", "SC"), (9, 5, "11", "DC"),
    (4, 6, "40", "OK"), (5, 6, "22", "LA"), (6, 6, "28", "MS"), (7, 6, "01", "AL"),
    (8, 6, "13", "GA"),
    (0, 7, "15", "HI"), (4, 7, "48", "TX"), (8, 7, "12", "FL"),
];
