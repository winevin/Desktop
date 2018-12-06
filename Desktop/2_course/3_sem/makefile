CFLAGS=-g -fsanitize=address -fsanitize=undefined -O2 -std=gnu11 -Wall -Wno-pointer-sign -Wformat -Wformat-overflow -Wformat-security -Wnull-dereference -Wignored-qualifiers -Wshift-negative-value -Wswitch-default -Wduplicated-branches -Wduplicated-branches -Wfloat-equal -Wshadow -Wpointer-arith -Wpointer-compare -Wtype-limits -Wwrite-strings -Wdangling-else -Wempty-body -Wlogical-op -Wstrict-prototypes -Wold-style-declaration -Wold-style-definition -Wmissing-parameter-type -Wmissing-field-initializers -Wnested-externs -Wvla-larger-than=4096 -lm -Werror=pointer-sign -Werror=format -Werror=unused-variable
all:
	clear
	gcc $(CFLAGS) $(NAME) -o main
clean:
	rm -rf main*
wifi:
	sudo rfkill block wifi
	sudo rfkill unblock wifi
