C --------------------------------------------------------------------
C This file crates matlabdll.dll whih is ment to read the .m file
C
C Copyright (C) 2010-2014 European Commission
C
C This file is part of Program DMM
C
C DMM is free software developed at the Joint Research Centre of the
C European Commission: you can redistribute it and/or modify it under
C the terms of the GNU General Public License as published by
C the Free Software Foundation, either version 3 of the License, or
C (at your option) any later version.
C
C DMM is distributed in the hope that it will be useful,
C but WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C GNU General Public License for more details.
C
C You should have received a copy of the GNU General Public License
C along with DMM.  If not, see <http://www.gnu.org/licenses/>.
C ---------------------------------------------------------------------
C -----------------------------------------
C To make dynamic the name of the .m file
C -----------------------------------------
	SUBROUTINE SETFILEM(string1,string2)
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:'setfilem_' :: SETFILEM
	CHARACTER*200 string1,string2,mfile,pathmfile
      COMMON /M/ mfile, pathmfile
      mfile     = string1
      pathmfile = string2
	RETURN
      END
