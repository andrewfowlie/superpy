TAGS = SM Hp Cha Neu Glu

all: f/bsll.h

f/bsll.h: f/%.h: $(foreach TAG,$(TAGS),m/%-$(TAG).F)
	mkdir -p f
	cat $^ > $@

m/bsll-%.F: bqll
	./bsll $(@:m/bsll-%.F=%)

clean:
	rm -fr m f

.PRECIOUS: %.F

