colorscheme ron
set hlsearch
set ts=4 sw=4 sts=4 et ai
filetype plugin on

"- Colorize OpenFOAM case files
execute pathogen#infect()
let g:foam256_use_custom_colors=1
set t_Co=256

"- Highlight long lines
autocmd FileType python,c,cpp,sh set colorcolumn=100

"- Turn off autoindent
autocmd FileType * setlocal formatoptions-=c formatoptions-=r formatoptions-=o

"- Highlight trailing whitespaces
highlight ExtraWhitespace ctermbg=red guibg=red
match ExtraWhitespace /\s\+$/
autocmd BufWinEnter * match ExtraWhitespace /\s\+$/
autocmd InsertEnter * match ExtraWhitespace /\s\+\%#\@<!$/
autocmd InsertLeave * match ExtraWhitespace /\s\+$/
autocmd BufWinLeave * call clearmatches()

"- Enables transparent pasting for vim < 8.0
"- Imported from https://stackoverflow.com/a/7053522
if &term =~ "xterm.*"
    let &t_ti = &t_ti . "\e[?2004h"
    let &t_te = "\e[?2004l" . &t_te
    function! XTermPasteBegin(ret)
        set pastetoggle=<Esc>[201~
        set paste
        return a:ret
    endfunction
    map <expr> <Esc>[200~ XTermPasteBegin("i")
    imap <expr> <Esc>[200~ XTermPasteBegin("")
    vmap <expr> <Esc>[200~ XTermPasteBegin("c")
    cmap <Esc>[200~ <nop>
    cmap <Esc>[201~ <nop>
endif
