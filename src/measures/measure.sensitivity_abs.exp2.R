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

source("src/measures/split-half.R")
source("src/measures/sensitivity_abs.R")

measure.sensitivity_abs.exp2.estimate <- function(X, n_t_)
{
	return( sensitivity_abs.estimate(X, n_t_, extrapolate.exp2) )
}
measure.sensitivity_abs.exp2.actual <- sensitivity_abs.actual
