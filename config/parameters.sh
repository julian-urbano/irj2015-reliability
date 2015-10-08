# Copyright (C) 2015  Julián Urbano <urbano.julian@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see http://www.gnu.org/licenses/.

# Create array with names of test collections in 'data/*.csv'.
# If you want to use specific files, create an array here:
# COLLECTIONS=(mycol1 mycol2)
COLLECTIONS=(data/*.csv)

# Clean collection names (remove path and file extension)
for I in ${!COLLECTIONS[*]}; do
	COLLECTION=${COLLECTIONS[$I]}
	COLLECTION=${COLLECTION##*/} # remove path
	COLLECTION=${COLLECTION%.*} # remove '.csv' extension
	COLLECTIONS[$I]=$COLLECTION
done

# Create array with names of measures in 'src/measures/measure.*.R'.
# If you want to use specific measures, create an array here:
# MEASURES=(mymeasure1 mymeasure2)
MEASURES=(src/measures/measure.*.R)

# Clean measure names (remove path and file extension)
for I in ${!MEASURES[*]}; do
        MEASURE=${MEASURES[$I]}
        MEASURE=${MEASURE##*/measure.} # remove path and 'measure.'
        MEASURE=${MEASURE%.*} # remove '.R'
        MEASURES[$I]=$MEASURE
done

# Create arrays with assumption values
NORMAL=(T F)
HOMOSCEDASTIC=(T F)
UNCORRELATED=(T F)
RANDOM_SAMPLING=(T F)
