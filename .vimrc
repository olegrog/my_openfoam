colorscheme ron
set hlsearch
set ts=4 sw=4 sts=4 et ai
filetype plugin indent on

"- Colorize OpenFOAM case files
execute pathogen#infect()
let g:foam256_use_custom_colors=1
set t_Co=256

"- Highlight long lines
autocmd FileType python,c,cpp,sh set colorcolumn=100

"- Highlight trailing whitespaces
highlight ExtraWhitespace ctermbg=red guibg=red
match ExtraWhitespace /\s\+$/
autocmd BufWinEnter * match ExtraWhitespace /\s\+$/
autocmd InsertEnter * match ExtraWhitespace /\s\+\%#\@<!$/
autocmd InsertLeave * match ExtraWhitespace /\s\+$/
autocmd BufWinLeave * call clearmatches()
