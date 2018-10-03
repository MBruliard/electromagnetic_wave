
tester:
	python3 test_matrices.py
	python3 test_projectors.py
	python3 test_energy.py

example:
	python3 test_eulerexpl.py
	python3 test_eulerimpl.py
	python3 test_rk2.py
	python3 test_rk4.py	
	
eulerexpl:
	python3 test_eulerexpl.py
	
eulerimpl:
	python3 test_eulerimpl.py

rk2:
	python3 test_rk2.py
	
rk4:
	python3 test_rk4.py


images_scheme:
	python3 comp_eulerexpl_analytic.py
	python3 comp_eulerimpl_analytic.py
	python3 comp_rk2_analytic.py
	python3 comp_rk4_analytic.py

order:
	python3 order_euler_expl.py
	python3 order_rk2.py

energy:
	python3 test_energy.py


# ------------------- OPTIONS ------------------------- #
clean:
	rm -Rf __pycache__
	
	
images_clean:
	rm Images/rk2/*
	rm Images/rk4/*
	rm Images/implicit/*
	rm Images/explicit/*.png
