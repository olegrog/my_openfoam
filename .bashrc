# .bashrc is sourced in Ubuntu

# Some of bash colors
export RED='\033[0;31m'
export YELLOW='\033[0;33m'
export BBLUE='\033[1;34m'
export NC='\033[0m'

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
    result=$(locate "/$1" | grep $WM_PROJECT_DIR | grep -v lnInclude)
    if [[ $? > 0 ]]; then
        echo -e "${RED}Missing filename!${NC}"
    elif [[ "$(wc -w <<< "$result")" > 1 ]]; then
        echo -e "${RED}There are several files with the same name."
        echo -e "Specify the unique path using the parent directory.${NC}"
        echo "$result" | tr ' ' '\n'
    else
        less -R $result
    fi
}

# Aliases
alias g='git'
alias grep='grep --color=auto'
alias ls='ls --color=auto -G'
alias ll='ls -l'
alias diff='colordiff'
alias solv='cd $FOAM_SOLVERS'
alias src='cd $FOAM_SRC'
alias tut='cd $FOAM_TUTORIALS'
