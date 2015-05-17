all:
	cd src ; cmake . 
	$(MAKE) -C src
	$(MAKE) -C timetests/
