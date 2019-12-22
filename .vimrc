colorscheme ron
set ts=4 sw=4 sts=4 et ai
set hlsearch

"- Colorize OpenFOAM case files
execute pathogen#infect()
filetype plugin indent on
let g:foam256_use_custom_colors=1
set t_Co=256

"- Highlight trailing whitespaces
highlight ExtraWhitespace ctermbg=red guibg=red
match ExtraWhitespace /\s\+$/
autocmd BufWinEnter * match ExtraWhitespace /\s\+$/
autocmd InsertEnter * match ExtraWhitespace /\s\+\%#\@<!$/
autocmd InsertLeave * match ExtraWhitespace /\s\+$/
autocmd BufWinLeave * call clearmatches()
