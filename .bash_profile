# .bash_profile is sourced in Centos

_source_if_exists() {
    [ -f "$1" ] && source "$@"
}

_source_if_exists "${HOME}/.bashrc"
