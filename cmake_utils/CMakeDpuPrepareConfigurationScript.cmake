file(COPY ${source}/${script} DESTINATION ${destination})
file(APPEND ${destination}/${script} "\nsave ${destination}/${output}\nquit\n")