function stim = df_Check_And_Load(file_name)
%   stim = df_Check_And_Load(file_name)
%       -- Check file type of given file_name 
%       -- Load file and save to the particular variable stim
% 
%             STRFPAK: STRF Estimation Software
% Copyright �2003. The Regents of the University of California (Regents).
% All Rights Reserved.
% Created by Theunissen Lab and Gallant Lab, Department of Psychology, Un
% -iversity of California, Berkeley.
%
% Permission to use, copy, and modify this software and its documentation
% for educational, research, and not-for-profit purposes, without fee and
% without a signed licensing agreement, is hereby granted, provided that
% the above copyright notice, this paragraph and the following two paragr
% -aphs appear in all copies and modifications. Contact The Office of Tec
% -hnology Licensing, UC Berkeley, 2150 Shattuck Avenue, Suite 510,
% Berkeley, CA 94720-1620, (510) 643-7201, for commercial licensing
% opportunities.
%
%IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
%SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
%ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
%REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
%LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
%PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY,
%PROVIDED HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PRO
%-VIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

%Created by JXZ, Sept, 2002.

% Check if the file exists
if exist(file_name, 'file')==0
    filename = sprintf('The given file: %s does not exist.', file_name);
    errordlg(filename, 'File error', 'modal');
    return;
end

[path,name,ext] = fileparts(file_name);
switch ext
case {'.dat', '.txt'}
     stim = load (file_name);
case {'.mat'}
     stimMat = load(file_name);
     % Validate the MAT-file
     flds = fieldnames(stimMat);
     if (length(flds) == 1)
         stim = getfield(stimMat, char(flds{1}));
     elseif (length(flds) > 1)
         for ii = 1:length(flds)
             stim{ii} = getfield(stimMat, char(flds{ii}));
         end
     end
otherwise
     errordlg('Wrong data file type.', 'File Type Error', 'modal')
end

% ===========================================================
% END OF CHECK_AND_LOAD
% ===========================================================
