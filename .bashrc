# .bashrc is sourced in Ubuntu

# Some of bash colors
declare -xr RED='\033[0;31m'
declare -xr YELLOW='\033[0;33m'
declare -xr BBLUE='\033[1;34m'
declare -xr NC='\033[0m'

# Command line
_parse_error_code() {
    local code=$?
    [[ $code > 0 ]] || return
    echo "[$code] "
}
# \[..color..\] should be used to preserve correct shell line width
export PS1="\[$RED\]\$(_parse_error_code)\[$YELLOW\]\h\[$NC\]:\[$BBLUE\]\w\[$NC\]\$ "

# Bash history
export HISTTIMEFORMAT="%d/%m/%y %T "
export HISTSIZE=10000
export HISTFILESIZE=10000
export HISTIGNORE="ls:ll:cd:cd -"
export HISTCONTROL="ignorespace:erasedups"
shopt -s histappend

# Colorize less output
export LESSOPEN="| $(which source-highlight) -i %s -f esc256 --style-file esc256.style"
export LESS=" -R"

# Search among OpenFOAM source
s() {
    less `locate "/$1" | grep $WM_PROJECT_DIR | grep -v applications | head -1`
}

# Aliases
alias g='git'
alias grep='grep --color=auto'
alias ls='ls --color=auto -G'
alias ll='ls -l'
alias diff='colordiff'
