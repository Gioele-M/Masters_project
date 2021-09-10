import threading
import time
import queue

def do_something(x,y):
    return f'you sent {x} and {y}'


if __name__ == '__main__':
    
    start = time.perf_counter()
    list_of_stuff = [4,4,3,2,1,2,4,2,4,24,2,42]
    other_thing = 'text'
    threads = []
    results = []

    for thing in list_of_stuff:
        t = threading.Thread(target=lambda: results.append(do_something(thing, other_thing)))
        t.start()
        threads.append(t)
    
    for thread in threads:
        thread.join()

    print(results)
    