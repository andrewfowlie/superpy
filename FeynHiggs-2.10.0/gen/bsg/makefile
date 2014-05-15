TAGS = SM Hp Cha Neu Glu

all: f/bsgamma.h f/bsglu.h

f/bsgamma.h f/bsglu.h: f/%.h: $(foreach TAG,$(TAGS),m/%-$(TAG).F)
	mkdir -p f
	cat $^ > $@

m/bsgamma-%.F: bsg
	./bsgamma $(@:m/bsgamma-%.F=%)

m/bsglu-%.F: bsg
	./bsglu $(@:m/bsglu-%.F=%)

clean:
	rm -fr m f

.PRECIOUS: %.F

