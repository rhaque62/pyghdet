""" A python file to test the pytorn package."""
from pyghdet import *

unq_species = ['a', 'b', 'c', 'd']
sus_species = ['e']

def test_spcom1():
    assert spcomb(unq_species , sus_species) == "Error:The provided suspected hybrid/s ['e'] is/are not in the list of species in the data!"

def test_spcom2():
    assert len(spcomb(unq_species, unq_species)) == 12

def test_cct1():
    assert cct([0, 0.05, 0.01, 0.99]) == 0

def test_cct2():
    assert cct([0.04, 0.001, 1]) == 1

def test_cct3():
    assert cct([0.01, 0.05, 0.99, 0.001]) <= 1
    assert cct([0.01, 0.05, 0.99, 0.001]) >= 0

def test_cct4():
    assert cct([0.01, 0.05, 'a']) == "Warning: The individual tests produced p-values containing non-numeric character! Failed to test the global null hypothesis"

def test_cct5():
    assert cct([0.01, -0.05, 0.09]) == "Warning: All the individual p-values must be between 0 and 1! Failed to test the global null hypothesis"
    assert cct([0.01, 1.2, 0.09]) == "Warning: All the individual p-values must be between 0 and 1! Failed to test the global null hypothesis"
    assert cct([0.01, -0.05, 0.09, 1.5]) == "Warning: All the individual p-values must be between 0 and 1! Failed to test the global null hypothesis"

def test_cct6():
    assert cct([0.01, 0, 0.99, 1]) == "Error: cannot have both 0 and 1 p-values!"


def test_cct7():
    assert cct([0.01, 0.05, 0.99], [1,2]) == "Error: weights and pvlaues should be same length!"
    assert cct([0.01, 0.05, 0.99], [1,2,3,4]) == "Error: weights and pvlaues should be same length!"


def test_cct8():
    assert cct([0.01, 0.05, 0.99], [-1, 2, 3]) == "Error: All the weights must be positive!"


def test_mcm():
    assert mcm([0.01, 0.05, 0.99, 0.001]) <= 1
    assert mcm([0.01, 0.05, 0.99, 0.001]) >= 0

def test_cmc():
    assert cmc([0.01, 0.05, 0.99, 0.001]) <= 1
    assert cmc([0.01, 0.05, 0.99, 0.001]) >= 0

def test_indiv():
    res = comb_indiv("data.txt", "map.txt", "out", 16, 4, 50000, alpha = 0.05)
    assert 0<= res.p_value <= 1

    res2 =  comb_indiv("data.txt", "map.txt", "out", 16, 4, 50000, alpha = 0.10)
    assert 0<= res2.p_value <= 1
    assert res.p_value == res2.p_value

    res3 = comb_indiv("data.txt", "map.txt", "out", 16, 4, 50000, alpha = 0.10, remove_amb_site = True)
    assert 0<= res3.p_value <= 1


def test_species():
    res = comb_species("data.txt", "map.txt", "out", 16, 4, 50000, alpha = 0.05)
    assert 0<= res.p_value <= 1

    res2 =  comb_species("data.txt", "map.txt", "out", 16, 4, 50000, alpha = 0.10)
    assert 0<= res2.p_value <= 1
    assert res.p_value == res2.p_value

    res3 = comb_species("data.txt", "map.txt", "out", 16, 4, 50000, alpha = 0.10, remove_amb_site = True)
    assert 0<= res3.p_value <= 1


def test_indiv2():
    res = comb_indiv("data.txt", "map.txt", "out", 16, 4, 50000,sus_hyb=['sp1'], alpha = 0.05)
    assert 0<= res.p_value <= 1

    res2 =  comb_indiv("data.txt", "map.txt", "out", 16, 4, 50000,sus_hyb=['sp1'], alpha = 0.10)
    assert 0<= res2.p_value <= 1
    assert res.p_value == res2.p_value

    res3 = comb_indiv("data.txt", "map.txt", "out", 16, 4, 50000,sus_hyb=['sp1'], alpha = 0.10, remove_amb_site = True)
    assert 0<= res3.p_value <= 1


def test_species2():
    res = comb_species("data.txt", "map.txt", "out", 16, 4, 50000,sus_hyb=['sp1'], alpha = 0.05)
    assert 0<= res.p_value <= 1

    res2 =  comb_species("data.txt", "map.txt", "out", 16, 4, 50000,sus_hyb=['sp1'], alpha = 0.10)
    assert 0<= res2.p_value <= 1
    assert res.p_value == res2.p_value

    res3 = comb_species("data.txt", "map.txt", "out", 16, 4, 50000,sus_hyb=['sp1'], alpha = 0.10, remove_amb_site = True)
    assert 0<= res3.p_value <= 1


def test_indiv3():
    res = comb_indiv("data.txt", "map.txt", "out", 16, 4, 50000,sus_hyb=['sp8'], alpha = 0.05)
    assert res == "Error:The provided suspected hybrid/s ['sp8'] is/are not in the list of species in the data!"



def test_species3():
    res = comb_species("data.txt", "map.txt", "out", 16, 4, 50000,sus_hyb=['sp8'], alpha = 0.05)
    assert res == "Error:The provided suspected hybrid/s ['sp8'] is/are not in the list of species in the data!"
