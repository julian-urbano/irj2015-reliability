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

# Check that we are called with a command. If not, exit.
if [ -z "$COMMAND" ]; then
	echo "Don't execute this script directly. Check the documentation."
	exit 1
fi

source config/environment.sh
source config/parameters.sh

# Check whether we have to loop collections in this command.
# If not, set to "" to iterate at least once.
if [ "$LOOP_COLLECTIONS" != true ]; then
	COLLECTIONS=("")
fi

# Check whether we have to loop measures in this command.
# If not, set to "" to iterate at least once.
if [ "$LOOP_MEASURES" != true ]; then
	MEASURES=("")
fi

# Check whether we have to loop assumptions in this command.
# If not, set to "" to iterate at least once.
if [ "$LOOP_ASSUMPTIONS" != true ]; then
	NORMAL=("")
	HOMOSCEDASTIC=("")
	UNCORRELATED=("")
	RANDOM_SAMPLING=("")
fi

# Create dir for qsub's out and err files
mkdir -p scratch/qout

for MEASURE in "${MEASURES[@]}"; do
	for COLLECTION in "${COLLECTIONS[@]}"; do
		for N in "${NORMAL[@]}"; do
			for H in "${HOMOSCEDASTIC[@]}"; do
				for U in "${UNCORRELATED[@]}"; do
					for R in "${RANDOM_SAMPLING[@]}"; do
						if [ "$SGE" == true ]; then
							qsub -v COMMAND=$COMMAND -v MEASURE=$MEASURE -v COLLECTION=$COLLECTION \
								-v NORMAL=$N -v HOMOSCEDASTIC=$H -v UNCORRELATED=$U -v RANDOM_SAMPLING=$R \
								-N "$COLLECTION$MEASURE$N$H$U$R($COMMAND)" \
								bin/qsub.sub
							sleep $SGE_WAIT
						else
							echo "$RSCRIPT" src/${COMMAND}.R $MEASURE $COLLECTION $N $H $U $R
							"$RSCRIPT" src/${COMMAND}.R $MEASURE $COLLECTION $N $H $U $R
						fi
					done
				done
			done
		done
	done
done
