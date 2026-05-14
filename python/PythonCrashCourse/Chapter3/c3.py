#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Sat 24. Jan 2026
#
# Purpose:      Go throu 
# =============================================================================

friend_list =["Victor", "Anzo", "Enyel", "Gabo", "Egor", "Pol", "Ailis", "Lucy"]

def list_friends():
    print(f"\tTotal number of friends: {len(friend_list)}")
    print( "---------------------------------------------")
    for friend_num in range(len(friend_list)):
        print(f"\tFriend number {friend_num+1} is {friend_list[friend_num]}")
    print( "---------------------------------------------")

def main() -> None:
    list_friends()
    check_new_friend = input("\nDo you have a new friend?\n").strip().title()
    if check_new_friend == "Yes": 
        who_new_friend = input("\nWho is it, liar?\n").strip().title()
        friend_list.append(who_new_friend)
        list_friends()
    Victor_returned_book = input("\nHave vic returned the book?\n").strip().title()
    if Victor_returned_book =="No": 
        del friend_list[0]
        list_friends()
        Victor_returned_book = input("\nHave vic returned the book?\n").strip().title()
        if Victor_returned_book=="Yes": 
            friend_list.insert(0,'Victor')
            list_friends()
        else:
            print("Que le jodan")
    else: 
        print("Bien por él")
    last_friend_made= friend_list.pop()
    print(f"\nThe last friend made was {last_friend_made}")
    check_removed_friend = input("\nDo you have a no longer friend?\n").strip().title()
    if check_removed_friend == "Yes": 
        who_nolonger_friend = input("\nWho is it, king?\n").strip().title()
        friend_list.remove(who_nolonger_friend)
        list_friends()

    print("\nThe real OGs:\n")
    friend_list.sort(reverse=True)
    for n in range(len(friend_list)):
        friend_list[n] = friend_list[n].upper() 
    list_friends()


if __name__ == "__main__":
    main()

